import gzip
import logging
import pathlib
import time
from collections import defaultdict
from datetime import datetime
from typing import Generator

import yaml
from ga4gh.vrs._internal.models import Allele  # noqa  F401 'Allele' private member
from tqdm import tqdm

import vrs_anvil
from vrs_anvil import Manifest, ThreadedTranslator, generate_gnomad_ids
from vrs_anvil.collector import collect_manifest_urls

_logger = logging.getLogger("vrs_anvil.annotator")


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
    """Return a gnomad expression generator for each line in the vcf."""
    all_done = False
    total_lines = 0
    for work_file in tqdm(
        _work_file_generator(manifest), total=len(manifest.vcf_files)
    ):
        if all_done:
            break
        line_number = 0
        if "gz" in str(work_file):
            f = gzip.open(work_file, "rt")
        else:
            f = open(work_file, "r")
        with f:
            key = str(work_file)
            metrics[key]["status"] = "started"
            metrics[key]["start_time"] = time.time()
            metrics[key]["successes"] = 0
            metrics[key]["metakb_hits"] = 0

            for line in f:
                if line.startswith("#"):
                    continue
                for gnomad_id in generate_gnomad_ids(
                    line, compute_for_ref=manifest.compute_for_ref
                ):
                    yield {"fmt": "gnomad", "var": gnomad_id}, work_file, line_number

                line_number += 1
                total_lines += 1

                if manifest.limit and line_number > manifest.limit:
                    _logger.info(f"Limit of {manifest.limit} reached, stopping")
                    all_done = True
                    break

            _logger.info(f"Setting metrics for {work_file}")
            metrics[key]["status"] = "finished"
            metrics[key]["end_time"] = time.time()
            metrics[key]["line_count"] = line_number
            metrics[key]["elapsed_time"] = (
                metrics[key]["end_time"] - metrics[key]["start_time"]
            )
            # TODO - should we delete the file (if its not a symlink) after we are done with it?

    _logger.info(
        f"_vcf_generator: Finished processing all files in the manifest {total_lines} lines processed."
    )


def _vrs_generator(manifest: Manifest) -> Generator[dict, None, None]:
    """Return a generator for the VRS ids."""
    tlr = ThreadedTranslator(normalize=manifest.normalize)
    c = 0
    for result in tlr.threaded_translate_from(
        generator=tqdm(_vcf_generator(manifest), total=manifest.estimated_vcf_lines),
        num_threads=manifest.num_threads,
    ):
        yield result
        c += 1
    _logger.info(
        f"_vrs_generator: Finished processing all vrs results in the manifest {c} results processed."
    )


def vrs_ids(allele: Allele) -> list[str]:
    """Return a list of VRS ids from an allele."""
    return [allele.id]  # , allele.location.id, allele.location.sequence_id]


def annotate_all(manifest: Manifest, max_errors: int) -> pathlib.Path:
    """Annotate all the files in the manifest. Return a file with metrics."""

    # set the manifest in a well known place, TODO: is this really necessary
    _logger.info("annotate_all: Starting.")
    vrs_anvil.manifest = manifest
    metakb_proxy = vrs_anvil.MetaKBProxy(
        metakb_path=pathlib.Path(manifest.metakb_directory)
    )
    _logger.info("annotate_all: completed metakb init.")

    metrics["total"]["start_time"] = time.time()
    total_errors = 0
    for result_dict in _vrs_generator(manifest):
        assert result_dict is not None, "result_dict is None"
        assert isinstance(result_dict, dict), "result_dict is not a dict"
        for k in ["file", "line"]:
            assert (
                k in result_dict
            ), f"metrics tracking from caller {k} not in result_dict"
            assert (
                result_dict[k] is not None
            ), f"metrics tracking from caller {k} is None"

        key = str(result_dict["file"])
        if "error" in result_dict:
            errors = metrics[key]["errors"]
            if result_dict["error"] not in errors:
                errors[result_dict["error"]] = 0
            errors[result_dict["error"]] += 1
            total_errors += 1
            if total_errors > max_errors:
                break
        else:
            result = result_dict.get("result", None)
            assert isinstance(
                result, Allele
            ), f"result is not the expected Pydantic Model {type(result)} {result_dict.keys()}"
            metrics[key]["successes"] += 1
            # check metaKB cache, TODO - it would be nice if we had the metakb.study.id and added that to result_dict
            if any([metakb_proxy.get(_) for _ in vrs_ids(result)]):
                _logger.info(f"VRS id {result.id} found in metakb. {result_dict}")
                metrics[key]["metakb_hits"] += 1

    _logger.info("annotate_all: Finished processing results.")

    metrics["total"]["end_time"] = time.time()
    metrics["total"]["elapsed_time"] = (
        metrics["total"]["end_time"] - metrics["total"]["start_time"]
    )
    metrics["total"]["successes"] = sum(
        [metrics[key].get("successes", 0) for key in metrics.keys() if key != "total"]
    )
    metrics["total"]["errors"] = sum(
        [
            sum(metrics[key]["errors"].values())
            for key in metrics.keys()
            if key != "total"
        ]
    )

    _logger.info("annotate_all: Finished calculating metrics.")

    # Append timestamp to filename
    timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    metrics_file = (
        pathlib.Path(manifest.state_directory) / f"metrics_{timestamp_str}.yaml"
    )
    with open(metrics_file, "w") as f:
        # clean up the recursive dict into a plain old dict so that it serialized to yaml neatly
        for k, v in metrics.items():
            metrics[k] = dict(v)
            if "errors" in metrics[k] and k != "total":
                metrics[k]["errors"] = dict(metrics[k]["errors"])
        yaml.dump(dict(metrics), f)

    _logger.info("annotate_all: Finished writing metrics.")

    return metrics_file
