def validate_threaded_result(result, errors, validate_passthrough=False):
    """Test helper, a utility method to validate threaded lookup result."""
    assert result is not None, "result_dict is None"
    from vrs_anvil.translator import VCFItem
    assert isinstance(result, VCFItem), "result_dict is not a dict"
    allele = result.result
    if isinstance(allele, dict):
        _ = allele  # handle dummy results
    else:
        _ = allele.model_dump(exclude_none=True)
    assert _ is not None, "result from allele lookup is missing"
    if "error" in _:
        errors.append(_)
        return
    for k in ["location", "state", "type"]:
        assert k in _, f"{k} not in result from allele lookup"
    for k in ["file_name", "line_number"]:
        assert getattr(result, k), f"metrics tracking from caller {k} not in result_dict {result}"
        if validate_passthrough:
            assert (
                    result[k] is not None
            ), f"metrics tracking from caller {k} is None"
