import pathlib

import pytest

from vrs_anvil import meta_kb_ids, MetaKBProxy


@pytest.fixture
def expected_vrs_id_count():
    return 2970


def test_metakb_ids(metakb_directory, expected_vrs_id_count, initialized_manifest):
    """Test metakb ids."""
    vrs_ids = [_ for _ in meta_kb_ids(metakb_directory)]
    vrs_count = len(vrs_ids)
    assert vrs_count >= expected_vrs_id_count, f"Not enough VRS ids found in metakb {vrs_count} {vrs_ids}"
    metakb_proxy = MetaKBProxy(metakb_path=pathlib.Path(metakb_directory))

    for _ in vrs_ids:
        assert metakb_proxy.get(_), f"VRS id {_} not found in cache {_}"

    hits, misses = metakb_proxy._cache.stats()
    assert misses == 0, f"Misses found in cache {misses}"
