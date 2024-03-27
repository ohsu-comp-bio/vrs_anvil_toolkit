import pytest

from vrs_anvil import METAKB_API, query_metakb


@pytest.fixture()
def vrs_id_1():
    return "ga4gh:VA.SOEVGpU16hxYQtJNeRyfq0V-B0rSOGK-"


def test_successful_query(vrs_id_1):
    response = query_metakb(vrs_id_1)
    assert response is not None, (
        f"unsuccessful query for {vrs_id_1}, "
        + f"ensure VRS ID digests have not changed at {METAKB_API}"
    )

    assert len(response["study_ids"]) != 0, (
        "no study ids found, "
        + f"ensure response format or VRS ID digests have not changed at {METAKB_API}"
    )
