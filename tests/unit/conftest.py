import pytest


@pytest.fixture
def seqrepo_dir():
    with open(".env") as f:
        for line in f:
            # Ignore comments and empty lines
            if line.strip() and not line.strip().startswith('#'):
                key, value = line.strip().split('=', 1)
                if key == "SEQREPO_ROOT":
                    return value + "/latest"
