from click.testing import CliRunner

from vrs_anvil import Manifest
from vrs_anvil.cli import cli


def test_cli_version():
    """Ensure version printed"""
    runner = CliRunner()
    result = runner.invoke(cli, "--version".split())
    expected_strings = ["version"]
    print(result.output)
    for expected_string in expected_strings:
        assert (
            expected_string in result.output
        ), f"Should have printed {expected_string}"


def test_loading_manifest(manifest_path):
    """Ensure manifest parsed and stored in context"""
    runner = CliRunner()
    print(manifest_path)
    result = runner.invoke(cli, f"--manifest {manifest_path} --verbose".split())
    expected_strings = Manifest.model_fields.keys()
    print(result.output)
    for expected_string in expected_strings:
        assert (
            expected_string in result.output
        ), f"Should have printed {expected_string}"
