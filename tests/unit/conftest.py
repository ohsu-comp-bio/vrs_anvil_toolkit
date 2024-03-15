import logging
import pathlib

import pytest
import yaml

import vrs_anvil
from vrs_anvil import caching_allele_translator_factory, ThreadedTranslator, Manifest

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


@pytest.fixture
def testing_manifest(manifest_path: pathlib.Path, tmp_path) -> Manifest:
    """Set the testing manifest for this run."""
    with open(manifest_path, "r") as stream:
        manifest_dict = yaml.safe_load(stream)
        for _ in ["work_directory", "cache_directory", "state_directory"]:
            manifest_dict[_] = str(tmp_path / manifest_dict[_])
        manifest = Manifest.parse_obj(manifest_dict)
    return manifest


@pytest.fixture
def testing_manifest_path(testing_manifest: Manifest, tmp_path) -> pathlib.Path:
    """Set the testing manifest for this run."""
    _ = tmp_path / "manifest.yml"
    with open(_, "w") as stream:
        yaml.dump(testing_manifest.model_dump(), stream)
    return _


@pytest.fixture
def metakb_directory(testing_manifest):
    """Where the CDM files are located."""
    return testing_manifest.metakb_directory


@pytest.fixture
def initialized_manifest(testing_manifest) -> bool:
    """state change, the manifest is initialized."""
    vrs_anvil.manifest = testing_manifest
    return True


@pytest.fixture
def seqrepo_dir(testing_manifest):
    """Return the seqrepo directory as fixture."""
    return testing_manifest.seqrepo_directory


@pytest.fixture
def my_translator(initialized_manifest):
    """Return a single translator instance."""
    # relies on testing_manifest setting the vrs_anvil.manifest
    return caching_allele_translator_factory()


@pytest.fixture
def threaded_translator(initialized_manifest):
    """Return a "thread manager", a pool of translator instances."""
    return ThreadedTranslator(normalize=False)


@pytest.fixture()
def num_threads(testing_manifest):
    """Return the number of threads to use for testing."""
    return vrs_anvil.manifest.num_threads


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
