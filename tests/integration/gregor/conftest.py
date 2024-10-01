import os
import pytest


@pytest.fixture
def vcf_path():
    """Return a remote VCF URL."""
    assert (
        "VCF_PATH" in os.environ
    ), "no VCF specified to run test on, please export VCF_PATH with the VCF to use for testing"
    return os.environ["VCF_PATH"]
