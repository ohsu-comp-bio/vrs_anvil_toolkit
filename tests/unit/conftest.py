import json
import logging
import pathlib
from typing import Generator

import pytest

from vrs_anvil import seqrepo_dir, caching_allele_translator_factory, ThreadedTranslator, find_items_with_key

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)


@pytest.fixture
def seqrepo_dir():
    """Return the seqrepo directory as fixture."""
    return seqrepo_dir()


@pytest.fixture
def my_translator():
    """Return a single translator instance."""
    return caching_allele_translator_factory()


@pytest.fixture
def threaded_translator():
    """Return a "thread manager", a pool of translator instances."""
    return ThreadedTranslator(normalize=False)


@pytest.fixture()
def num_threads():
    """Return the number of threads to use for testing."""
    return 20


@pytest.fixture
def thousand_genome_vcf() -> pathlib.Path:
    """Return a path to a 1000g vcf."""
    _ = pathlib.Path("tests/fixtures/1kGP.chr1.1000.vcf")
    assert _.exists()
    return _
