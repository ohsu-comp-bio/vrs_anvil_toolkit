import os
import pytest

from click.testing import CliRunner
from vrs_anvil import Manifest
from vrs_anvil.cli import cli


@pytest.fixture
def recent_timestamp():
    return "20240507_203113"


@pytest.fixture
def scatter_path():
    return "tests/fixtures/scatter/"


@pytest.fixture
def manifest_path_ps():
    return "tests/fixtures/manifest_ps.yaml"


@pytest.fixture
def num_vcfs(testing_manifest):
    return len(testing_manifest.vcf_files)


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


def test_loading_manifest(manifest_path):
    """Test that manifest can be parsed and stored in context"""
    runner = CliRunner()
    result = runner.invoke(cli, f"--manifest {manifest_path} --verbose")
    expected_strings = Manifest.model_fields.keys()
    print(result.output)
    for expected_string in expected_strings:
        assert (
            expected_string in result.output
        ), f"Should have printed {expected_string}"


# def test_annotate_scatter(manifest_path, testing_manifest):
#     """Test that manifest outputs the right files"""
#     # run the command
#     manifest = testing_manifest
#     runner = CliRunner()
#     with runner.isolated_filesystem() as test_dir:
#         shutil.copy(manifest_path, test_dir)
#         result = runner.invoke(cli, f"--manifest f{manifest_path}/{test_dir} annotate --scatter")

#     # assert metrics in

#     print(result.output)

#     # assert()


def test_ps_returns_recent_files(recent_timestamp, num_vcfs):
    """Test that vrs_anvil ps returns the most recent scatter command"""

    # make function call
    runner = CliRunner()
    os.chdir("tests/fixtures/ps")

    try:
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
    finally:
        os.chdir("../../..")
