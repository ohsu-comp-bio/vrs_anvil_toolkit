

def validate_threaded_result(result_dict, errors, validate_passthrough=False):
    """Test helper, a utility method to validate threaded lookup result."""
    assert result_dict is not None, "result_dict is None"
    assert isinstance(result_dict, dict), "result_dict is not a dict"
    result = result_dict.get('result', None)
    if isinstance(result, dict):
        _ = result  # handle dummy results
    else:
        _ = result.model_dump(exclude_none=True)
    assert _ is not None, "result from allele lookup is missing"
    if 'error' in _:
        errors.append(_)
        return
    for k in ['location', 'state', 'type']:
        assert k in _, f"{k} not in result from allele lookup"
    for k in ['file', 'line']:
        assert k in result_dict, f"metrics tracking from caller {k} not in result_dict"
        if validate_passthrough:
            assert result_dict[k] is not None, f"metrics tracking from caller {k} is None"
