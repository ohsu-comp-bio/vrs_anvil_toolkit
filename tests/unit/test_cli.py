from click.testing import CliRunner
from vrs_anvil.cli import cli


def test_version():
    """Ensure version printed"""
    runner = CliRunner()
    result = runner.invoke(cli, '--version'.split())
    expected_strings = ["version"]
    print(result.output)
    for expected_string in expected_strings:
        assert expected_string in result.output, f"Should have printed {expected_string}"


def test_manifest(manifest_path):
    """Ensure manifest parsed and stored in context"""
    runner = CliRunner()
    result = runner.invoke(cli, f'--manifest_path {manifest_path} --verbose'.split())
    expected_strings = ["work_directory='tmp/work'", "worker_count=20"]  # strings from manifest.yaml
    print(result.output)
    for expected_string in expected_strings:
        assert expected_string in result.output, f"Should have printed {expected_string}"
