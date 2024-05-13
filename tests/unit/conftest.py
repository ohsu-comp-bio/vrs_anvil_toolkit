import logging
import pathlib

import pytest
import yaml

from vrs_anvil.translator import Translator
from vrs_anvil import caching_allele_translator_factory, Manifest

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)


@pytest.fixture
def manifest_path() -> pathlib.Path:
    """Return a path to a manifest file."""
    _ = pathlib.Path("tests/fixtures/manifest.yaml").resolve()
    assert _.exists()
    return _


# @pytest.fixture(autouse=True)
# def change_test_dir(monkeypatch, tmp_path):
#     """Change to a temporary directory for testing."""
#     # see https://stackoverflow.com/questions/62044541/change-pytest-working-directory-to-test-case-directory
#     monkeypatch.chdir(tmp_path)


# note that tmp_path is a fixture provided by pytest
@pytest.fixture
def testing_manifest(manifest_path: pathlib.Path, tmp_path) -> Manifest:
    with open(manifest_path, "r") as stream:
        manifest_dict = yaml.safe_load(stream)

        # use pytest-provided relative path
        for dir_key in ["work_directory", "cache_directory", "state_directory"]:
            manifest_dict[dir_key] = str(tmp_path / manifest_dict[dir_key])

        manifest = Manifest.model_validate(manifest_dict)
    return manifest


@pytest.fixture
def metakb_directory(testing_manifest):
    """Where the CDM files are located."""
    return testing_manifest.metakb_directory


@pytest.fixture
def seqrepo_dir(testing_manifest):
    """Return the seqrepo directory as fixture."""
    return testing_manifest.seqrepo_directory


@pytest.fixture
def caching_translator(testing_manifest):
    """Return a single translator instance."""
    return caching_allele_translator_factory(
        seqrepo_directory=testing_manifest.seqrepo_directory
    )


@pytest.fixture
def translator():
    """Return a translator instance."""
    return Translator(normalize=False)


@pytest.fixture()
def num_threads(testing_manifest):
    """Return the number of threads to use for testing."""
    return testing_manifest.num_threads


@pytest.fixture
def thousand_genome_vcf() -> pathlib.Path:
    """Return a path to a 1000g vcf."""
    _ = pathlib.Path("tests/fixtures/1kGP.chr1.1000.vcf")
    assert _.exists()
    return _


@pytest.fixture
def python_source_directories() -> list[str]:
    """Directories to scan with flake8."""
    return ["vrs_anvil", "tests"]
