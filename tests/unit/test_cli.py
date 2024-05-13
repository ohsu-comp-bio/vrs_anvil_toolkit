import pytest

from click.testing import CliRunner
from vrs_anvil import Manifest
from vrs_anvil.cli import cli


@pytest.fixture
def recent_timestamp():
    return "20240507_203113"


@pytest.fixture
def manifest_path_ps():
    return "tests/fixtures/manifest_ps.yaml"


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


def test_ps_returns_recent_info(recent_timestamp, manifest_path_ps):
    """Test that vrs_anvil ps returns the most recent scatter command"""

    # make function call
    runner = CliRunner()
    result = runner.invoke(cli, f"--manifest {manifest_path_ps} ps")
    print(result.output)

    # check successful command and loaded most recent process file
    assert result.exit_code == 0, f"result failed with message: \n{result}"
    assert (
        "no scattered processes" not in result.output
    ), "no scattered processes located"
    assert recent_timestamp in result.output, "most recent date has not been chosen"

    # check metrics and manifest files located
    expected_num_processes = 0
    for i in range(expected_num_processes):
        assert (
            f"manifest_scattered_{recent_timestamp}_{i}.yaml" in result.output
        ), f"{i}th {expected_num_processes} manifest not found"

        assert (
            f"metrics_scattered_{recent_timestamp}_{i}.yaml" in result.output
        ), f"{i}th {expected_num_processes} metrics files not found"
