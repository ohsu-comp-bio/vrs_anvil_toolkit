import gzip
import logging
import pathlib
import time
from collections import defaultdict
from datetime import datetime
from typing import Generator

import yaml
from ga4gh.vrs import models as VRS
from tqdm import tqdm

import vrs_anvil
from vrs_anvil import Manifest, generate_gnomad_ids
from vrs_anvil.collector import collect_manifest_urls
from vrs_anvil.translator import Translator, VCFItem

_logger = logging.getLogger("vrs_anvil.annotator")

# enums for metrics
# TODO: do this for keys across files like "parameters" but also "fmt" and "line"
PARAMETERS = "parameters"
TOTAL = "total"
STATUS = "status"
SUCCESSES = "successes"
ERROR = "error"
METAKB_HITS = "metakb_hits"
MATCHES = "matches"
START_TIME = "start_time"
END_TIME = "end_time"
ELAPSED_TIME = "elapsed_time"
LINE_COUNT = "line_count"
VRS_OBJECT = "vrs_object"
TIMESTAMP = "timestamp_str"


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


def _vcf_item_generator(manifest: Manifest) -> Generator[tuple, None, None]:
    """Return a VCFItem for each line in the vcf."""
    total_lines = 0
    for work_file in tqdm(
        _work_file_generator(manifest),
        total=len(manifest.vcf_files),
        disable=manifest.disable_progress_bars,
    ):
        line_number = 0
        if "gz" in str(work_file):
            f = gzip.open(work_file, "rt")
        else:
            f = open(work_file, "r")
        with f:
            key = str(work_file)
            metrics[key][STATUS] = "started"
            metrics[key][START_TIME] = time.time()
            metrics[key][SUCCESSES] = 0
            metrics[key][METAKB_HITS] = 0

            for line in f:
                if line.startswith("#"):
                    continue

                line_number += 1
                total_lines += 1

                for gnomad_id in generate_gnomad_ids(
                    line, compute_for_ref=manifest.compute_for_ref
                ):
                    yield VCFItem(
                        fmt="gnomad",
                        var=gnomad_id,
                        file_name=work_file,
                        line_number=line_number,
                    )  # {"fmt": "gnomad", "var": gnomad_id}, work_file, line_number

                if manifest.limit and line_number > manifest.limit:
                    _logger.info(f"Limit of {manifest.limit} reached, stopping")
                    break

            _logger.info(f"Setting metrics for {work_file}")
            metrics[key][STATUS] = "finished"
            metrics[key][END_TIME] = time.time()
            metrics[key][LINE_COUNT] = line_number
            metrics[key][ELAPSED_TIME] = (
                metrics[key][END_TIME] - metrics[key][START_TIME]
            )

    _logger.info(
        f"_vcf_generator: Finished processing all files in the manifest {total_lines} lines processed."
    )


def _vrs_generator(manifest: Manifest) -> Generator[dict, None, None]:
    """Return a generator for the VRS ids."""
    tlr = Translator(normalize=manifest.normalize)
    c = 0
    for result in tlr.translate_from(
        generator=tqdm(
            _vcf_item_generator(manifest),
            total=manifest.estimated_vcf_lines,
            disable=manifest.disable_progress_bars,
        ),
        num_threads=manifest.num_threads,
    ):
        yield result
        c += 1
    _logger.info(
        f"_vrs_generator: Finished processing all vrs results in the manifest {c} results processed."
    )


def vrs_ids(allele: VRS.Allele) -> list[str]:
    """Return a list of VRS ids from an allele."""
    return [allele.id]  # , allele.location.id, allele.location.sequence_id]


def annotate_all(
    manifest: Manifest, max_errors: int, timestamp_str: str = None
) -> pathlib.Path:
    """Annotate all the files in the manifest. Return a file with metrics."""

    # set the manifest in a well known place, TODO: is this really necessary
    _logger.info("annotate_all: Starting.")
    vrs_anvil.manifest = manifest
    metakb_proxy = vrs_anvil.MetaKBProxy(
        metakb_path=pathlib.Path(manifest.metakb_directory),
        cache_path=pathlib.Path(manifest.cache_directory),
    )
    _logger.info("annotate_all: completed metakb init.")

    metrics[TOTAL][START_TIME] = time.time()
    total_errors = 0
    for result in _vrs_generator(manifest):
        assert result is not None, "result is None"
        assert isinstance(result, VCFItem), "result is not a VCFItem"

        file_path = str(result.file_name)

        if ERROR in result:
            errors = metrics[file_path][ERROR]
            if result[ERROR] not in errors:
                errors[result[ERROR]] = 0
            errors[result[ERROR]] += 1
            total_errors += 1
            if total_errors > max_errors:
                break
        else:
            allele_id = result.result

            metrics[file_path][SUCCESSES] += 1

            # check metaKB cache, TODO - it would be nice if we had the metakb.study.id and added that to result_dict
            if metakb_proxy.get(allele_id):
                _logger.info(f"VRS id {allele_id} found in metakb. {result}")

                # add vrs_id, allele_dict, actual evidence to this object as well (#3)
                metrics[file_path][MATCHES][allele_id] = {
                    "fmt": result.fmt,
                    "var": result.var,
                }

                metrics[file_path][METAKB_HITS] += 1

    _logger.info("annotate_all: Finished processing results.")

    metrics[TOTAL][TIMESTAMP] = timestamp_str
    metrics[TOTAL][END_TIME] = time.time()
    metrics[TOTAL][ELAPSED_TIME] = metrics[TOTAL][END_TIME] - metrics[TOTAL][START_TIME]
    metrics[TOTAL][SUCCESSES] = sum(
        [metrics[key].get(SUCCESSES, 0) for key in metrics.keys() if key != TOTAL]
    )
    metrics[TOTAL]["errors"] = sum(
        [sum(metrics[key]["errors"].values()) for key in metrics.keys() if key != TOTAL]
    )

    _logger.info("annotate_all: Finished calculating metrics.")

    # Append timestamp or suffix to filename
    if not timestamp_str:
        timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    metrics_file = (
        pathlib.Path(manifest.state_directory) / f"metrics_{timestamp_str}.yaml"
    )
    with open(metrics_file, "w") as f:
        # clean up the recursive dict into a plain old dict so that it serialized to yaml neatly
        for k, v in metrics.items():
            metrics[k] = dict(v)
            if k != TOTAL:
                if "errors" in metrics[k]:
                    metrics[k]["errors"] = dict(metrics[k]["errors"])
                if MATCHES in metrics[k]:
                    metrics[k][MATCHES] = dict(metrics[k][MATCHES])

        yaml.dump(dict(metrics), f)

    _logger.info("annotate_all: Finished writing metrics.")

    return metrics_file
