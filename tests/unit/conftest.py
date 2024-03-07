import logging
import pathlib

import pytest

from vrs_anvil import caching_allele_translator_factory, ThreadedTranslator

from vrs_anvil import seqrepo_dir as _seqrepo_dir

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)


@pytest.fixture
def seqrepo_dir():
    """Return the seqrepo directory as fixture."""
    return _seqrepo_dir()


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


@pytest.fixture
def manifest_path() -> pathlib.Path:
    """Return a path to a manifest file."""
    _ = pathlib.Path("tests/fixtures/manifest.yaml")
    assert _.exists()
    return _


@pytest.fixture
def python_source_directories() -> list[str]:
    """Directories to scan with flake8."""
    return ["vrs_anvil", "tests"]
