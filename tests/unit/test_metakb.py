import pathlib

import pytest

from vrs_anvil import metakb_ids, MetaKBProxy


@pytest.fixture
def expected_vrs_id_count():
    return 2986


@pytest.fixture
def expected_vrs_ids():
    return [
        "ga4gh:VA.SOEVGpU16hxYQtJNeRyfq0V-B0rSOGK-",
        "ga4gh:VA.fU8g-a8s0n-tFsj3-XbsuFe17MfySB4q",
        "ga4gh:VA.rRPCnh0XXjuePRGWerw6PhVXFYjhchwP",
    ]


def test_metakb_ids(
    metakb_directory, testing_manifest, expected_vrs_id_count, expected_vrs_ids
):
    """Test metakb ids."""
    vrs_ids = [vrs_id for vrs_id in metakb_ids(metakb_directory)]
    vrs_count = len(vrs_ids)
    assert (
        vrs_count >= expected_vrs_id_count
    ), f"Not enough VRS ids found in metakb {vrs_count} {vrs_ids}"

    metakb_proxy = MetaKBProxy(
        metakb_path=pathlib.Path(metakb_directory),
        cache_path=pathlib.Path(testing_manifest.cache_directory),
    )

    for id in vrs_ids:
        assert metakb_proxy.get(id), f"VRS id {id} not found in cache {id}"

    for id in expected_vrs_ids:
        assert metakb_proxy.get(id), f"Expected VRS id {id} not found in cache {id}"

    _, misses = metakb_proxy._cache.stats()
    assert misses == 0, f"Misses found in cache {misses}"
