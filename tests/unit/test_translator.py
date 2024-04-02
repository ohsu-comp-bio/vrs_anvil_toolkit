import logging
import pathlib
from typing import Generator

import pytest

from vrs_anvil.translator import threaded_translator, VCFItem

_logger = logging.getLogger("vrs_anvil.test_translator")


@pytest.fixture
def gnomad_csv() -> pathlib.Path:
    """Return a path to a gnomad csv file."""
    _ = pathlib.Path(
        "tests/fixtures/gnomAD_v4.0.0_ENSG00000012048_2024_03_04_18_33_26.csv"
    )
    assert _.exists()
    return _


def gnomad_ids(path, limit=None) -> Generator[VCFItem, None, None]:
    """Open the csv file and yield the first column 'gnomAD ID', skip the header"""
    c = 0
    with open(path, "r") as f:
        skip = True
        for line in f:
            if skip:
                skip = False
                continue
            _ = line.split(",")
            gnomad_id = _[0]
            # use allele number as the identifier
            allele_number = _[18]
            yield VCFItem("gnomad", gnomad_id, path, c, allele_number)
            c += 1
            if limit and c == limit:
                break


def test_threaded_translator(gnomad_csv):
    """Ensure threading works as expected."""

    limit = 2000
    num_worker_threads = 10

    generator_example = gnomad_ids(gnomad_csv, limit=limit)

    results = threaded_translator(generator_example, num_worker_threads)
    _logger.warning("test_threaded_translator: Starting test")

    # Add your assertions based on the expected results
    c = 0
    for _ in results:
        assert isinstance(_, VCFItem), "should get a VRS id"
        allele = _.result
        assert allele.id is not None, "allele.id is None"
        assert 'ga4gh:VA.' in allele.id, "allele.id is not a VRS id"
        c += 1

    if limit:
        assert c == limit, "did not get the expected number of results"
