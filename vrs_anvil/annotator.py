import pathlib
import time
from collections import defaultdict
from typing import Generator

import yaml
from ga4gh.vrs._internal.models import Allele  # noqa  F401 'Allele' private member
from tqdm import tqdm

import vrs_anvil
from vrs_anvil import Manifest, ThreadedTranslator, generate_gnomad_ids
from vrs_anvil.collector import collect_manifest_urls


def recursive_defaultdict():
    """Implicitly create an entry if a key is read that doesn’t yet exist, any level deep"""
    return defaultdict(recursive_defaultdict)


metrics = recursive_defaultdict()


def _work_file_generator(manifest: Manifest) -> Generator[pathlib.Path, None, None]:
    """Return a generator for the files in the manifest."""
    for work_file in collect_manifest_urls(manifest):
        work_file = pathlib.Path(work_file)
        assert work_file.exists(), f"File {work_file} does not exist"
        yield work_file


def _vcf_generator(manifest: Manifest) -> Generator[tuple, None, None]:
    """Return a generator lines in the vcf."""
    for work_file in tqdm(_work_file_generator(manifest)):
        line_number = 0
        with open(work_file, "r") as f:
            key = str(work_file)
            metrics[key]["status"] = 'started'
            metrics[key]["start_time"] = time.time()
            metrics[key]["successes"] = 0

            for line in f:
                line_number += 1
                if line.startswith("#"):
                    continue
                for gnomad_id in generate_gnomad_ids(line):
                    yield {"fmt": "gnomad", "var": gnomad_id}, work_file, line_number

            metrics[key]["status"] = 'finished'
            metrics[key]["end_time"] = time.time()
            metrics[key]["line_count"] = line_number
            metrics[key]["elapsed_time"] = metrics[key]["end_time"] - metrics[key]["start_time"]
            # TODO - should we delete the file (if its not a symlink) after we are done with it?


def _vrs_generator(manifest: Manifest) -> Generator[dict, None, None]:
    """Return a generator for the VRS ids."""
    tlr = ThreadedTranslator(normalize=manifest.normalize)
    for result in tlr.threaded_translate_from(generator=tqdm(_vcf_generator(manifest)),
                                              num_threads=manifest.num_threads):
        yield result


def annotate_all(manifest: Manifest, max_errors: int):
    """Annotate all the files in the manifest."""

    # set the manifest in a well known place, TODO: is this really necessary
    vrs_anvil.manifest = manifest

    metrics["total"]["start_time"] = time.time()
    total_errors = 0
    for result_dict in _vrs_generator(manifest):
        assert result_dict is not None, "result_dict is None"
        assert isinstance(result_dict, dict), "result_dict is not a dict"
        for k in ['file', 'line']:
            assert k in result_dict, f"metrics tracking from caller {k} not in result_dict"
            assert result_dict[k] is not None, f"metrics tracking from caller {k} is None"

        key = str(result_dict['file'])
        if 'error' in result_dict:
            errors = metrics[key]["errors"]
            if result_dict['error'] not in errors:
                errors[result_dict['error']] = 0
            errors[result_dict['error']] += 1
            total_errors += 1
            if total_errors > max_errors:
                break
        else:
            result = result_dict.get('result', None)
            assert isinstance(result, Allele), f"result is not the expected Pydantic Model {type(result)} {result_dict.keys()}"
            metrics[key]["successes"] += 1

    metrics["total"]["end_time"] = time.time()
    metrics["total"]["elapsed_time"] = metrics["total"]["end_time"] - metrics["total"]["start_time"]
    metrics["total"]["successes"] = sum([metrics[key].get("successes", 0) for key in metrics.keys() if key != "total"])
    metrics["total"]["errors"] = sum([sum(metrics[key]["errors"].values()) for key in metrics.keys() if key != "total"])

    metrics_file = pathlib.Path(manifest.state_directory) / "metrics.yaml"
    with open(metrics_file, "w") as f:
        # clean up the recursive dict into a plain old dict so that it serialized to yaml neatly
        for k, v in metrics.items():
            metrics[k] = dict(v)
            if 'errors' in metrics[k] and k != 'total':
                metrics[k]['errors'] = dict(metrics[k]['errors'])
        yaml.dump(dict(metrics), f)