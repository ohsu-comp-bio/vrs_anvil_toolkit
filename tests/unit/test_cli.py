import os
from pathlib import Path
import shutil
import pytest

from click.testing import CliRunner
from glob import glob
from unittest.mock import MagicMock, patch
from vrs_anvil import Manifest
from vrs_anvil.cli import cli

############
# FIXTURES #
############


@pytest.fixture
def num_vcfs(testing_manifest):
    return len(testing_manifest.vcf_files)


@pytest.fixture
def recent_timestamp():
    return "20240507_203113"


@pytest.fixture
def ps_dir():
    return "tests/fixtures/ps"


@pytest.fixture
def ps_manifest_path(ps_dir):
    return ps_dir + "/manifest.yaml"


@pytest.fixture
def mock_cli_manifest(tmp_path, monkeypatch, manifest_path, testing_manifest):
    '''Move relevant files and working dir into a temporary dir before testing cli commands'''

    # copy manifest and vcf files to tmp_path
    shutil.copy(manifest_path, tmp_path)
    print(tmp_path)

    for file_path in testing_manifest.vcf_files:
        path_arr = str(Path(file_path)).split("/")
        if len(path_arr) > 1:
            file_dir = "/".join(path_arr[:-1])
            os.makedirs(f"{tmp_path}/{file_dir}", exist_ok=True)
        shutil.copy(file_path, f"{tmp_path}/{file_path}")

    # create metakb cdm path
    metakb_rel_path = str(Path(testing_manifest.metakb_directory))
    metakb_test_dir = f"{tmp_path}/{str(metakb_rel_path)}"
    os.makedirs(metakb_test_dir)
    shutil.copytree(metakb_rel_path, metakb_test_dir, dirs_exist_ok=True)

    # change to temp testing directory
    monkeypatch.chdir(tmp_path)
    return testing_manifest


#########
# TESTS #
#########


@pytest.fixture
def suffix():
    return "arbitrary_suffix"


def test_cli_version():
    """Test that version can be called and printed"""
    runner = CliRunner()
    result = runner.invoke(cli, "--version")
    expected_strings = ["version"]
    print(result.output)
    for expected_string in expected_strings:
        assert (
            expected_string in result.output
        ), f"Should have printed {expected_string}"


def test_loading_manifest(mock_cli_manifest):
    # mock_cli_manifest used to set up directories

    """Test that manifest can be parsed and stored in context"""
    runner = CliRunner()
    result = runner.invoke(cli, "--verbose")
    expected_strings = Manifest.model_fields.keys()
    print(result.output)
    for expected_string in expected_strings:
        assert (
            expected_string in result.output
        ), f"Should have printed {expected_string}"


def test_using_suffix(mock_cli_manifest, suffix):
    # mock_cli_manifest used to set up directories
    """Test that the suffix is successfully added to the manifest string"""

    # stub the annotate_all function
    # note the path is to the namespace (module) being tested not the namespace its imported from
    mock = MagicMock()
    mock.return_value = "mocked_metrics.yaml"
    patcher = patch("vrs_anvil.cli.annotate_all", mock)
    patcher.start()

    runner = CliRunner()
    result = runner.invoke(cli, f"--suffix {suffix} annotate")

    print(result.output)
    assert suffix in result.output, f"Should have printed {suffix}"

    patcher.stop()


def test_annotate_scatter_outputs(mock_cli_manifest):
    """Test that manifest outputs the right files"""

    # setup
    manifest = mock_cli_manifest
    print("manifest", manifest)
    runner = CliRunner()

    # stub the annotate_all function
    mock = MagicMock()

    # run annotate scatter cmd
    with patch("vrs_anvil.cli.run_command_in_background", mock):
        result = runner.invoke(cli, "--verbose annotate --scatter")
    print(result.output)

    # make sure input files (manifest and scattered_processes) are created
    work_dir = str(Path(manifest.work_directory))
    expected_num_files = len(manifest.vcf_files)

    num_manifests = len(glob(f"{work_dir}/*manifest*.yaml"))
    assert (
        num_manifests == expected_num_files
    ), f"expected {expected_num_files} manifests, got {num_manifests}"

    process_files = len(glob(f"{work_dir}/*scattered_process*.yaml"))
    assert process_files == 1, f"expected 1 scattered process file, got {process_files}"

    # make sure main output log is created
    state_dir = str(Path(manifest.state_directory))
    num_log_files = len(glob(f"{state_dir}/*.log"))
    assert (num_log_files == 1), f"expected 1 log files, got {num_log_files}"


def test_ps_returns_recent_files(ps_dir, monkeypatch, recent_timestamp, num_vcfs):
    """Test that vrs_anvil ps returns the most recent scatter command"""

    # make function call
    runner = CliRunner()
    monkeypatch.chdir(ps_dir)

    result = runner.invoke(cli, "--manifest manifest.yaml ps")
    print(result.output)

    # check successful command and loaded most recent process file
    assert result.exit_code == 0, f"result failed with message: \n{result}"
    assert (
        "no scattered processes" not in result.output
    ), "no scattered processes located"
    assert recent_timestamp in result.output, "most recent date has not been chosen"

    # check metrics and manifest files located
    for i in range(num_vcfs):
        assert (
            f"manifest_scattered_{recent_timestamp}_{i}.yaml" in result.output
        ), f"manifest #{i} of {num_vcfs} not found"

        assert (
            f"metrics_scattered_{recent_timestamp}_{i}.yaml" in result.output
        ), f"metrics file #{i} of {num_vcfs} not found"
