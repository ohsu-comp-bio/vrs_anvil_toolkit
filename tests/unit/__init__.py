def validate_threaded_result(result, validate_passthrough=False):
    """Test helper, a utility method to validate threaded lookup result."""
    assert result is not None, "result_dict is None"
    from vrs_anvil.translator import VCFItem

    assert isinstance(result, VCFItem), "result_dict is not a dict"

    allele_id = result.result
    assert isinstance(
        allele_id, str
    ), f"translated VRS Allele ID is a {type(allele_id)} not a string"
    prefix, hash = allele_id.split(".")
    assert (
        prefix == "ga4gh:VA" and len(hash) == 32
    ), "VRS Allele ID format has changed, consult https://vrs.ga4gh.org/en/stable/impl-guide/computed_identifiers.html#identify"

    for k in ["file_name", "line_number"]:
        if validate_passthrough:
            assert (
                getattr(result, k) is not None
            ), f"metrics tracking from caller {k} is None"
