import pathlib
import time
from typing import Generator

import pytest

from tests.unit import validate_threaded_result
from vrs_anvil.translator import VCFItem

# see https://github.com/ga4gh/vrs-python/blob/main/tests/extras/test_allele_translator.py#L17

snv_inputs = {
    "hgvs": "NC_000019.10:g.44908822C>T",
    "beacon": "19 : 44908822 C > T",
    "spdi": "NC_000019.10:44908821:1:T",
    "gnomad": "19-44908822-C-T",
}

snv_output = {
    "id": "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt",
    "type": "Allele",
    "digest": "0AePZIWZUNsUlQTamyLrjm2HWUw2opLt",
    "location": {
        "id": "ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz",
        "type": "SequenceLocation",
        "digest": "wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
        },
        "start": 44908821,
        "end": 44908822,
    },
    "state": {"type": "LiteralSequenceExpression", "sequence": "T"},
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1373966/?new_evidence=true
deletion_inputs = {
    "hgvs": "NC_000013.11:g.20003097del",
    "spdi": ["NC_000013.11:20003096:C:", "NC_000013.11:20003096:1:"],
    "gnomad": "13-20003096-AC-A",
}

gnomad_deletion_output = {
    "id": "ga4gh:VA.agNTkBOVZ0HQP5iWJ2GoRovSaHUCyZHN",
    "type": "Allele",
    "digest": "agNTkBOVZ0HQP5iWJ2GoRovSaHUCyZHN",
    "location": {
        "id": "ga4gh:SL.30MgiQJ3IB1mk1wq4RyhUADF_V0vw0fD",
        "type": "SequenceLocation",
        "digest": "30MgiQJ3IB1mk1wq4RyhUADF_V0vw0fD",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
        },
        "start": 20003096,
        "end": 20003097,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        "length": 0,
        "sequence": "",
        "repeatSubunitLength": 1,
    },
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1687427/?new_evidence=true
insertion_inputs = {
    "hgvs": "NC_000013.11:g.20003010_20003011insG",
    "spdi": ["NC_000013.11:20003010::G", "NC_000013.11:20003010:0:G"],
    "gnomad": "13-20003010-A-AG",
}

gnomad_insertion_output = {
    "digest": "YSFR_-q58UM-bsgyq9ScHa4hfNOWGelM",
    "id": "ga4gh:VA.YSFR_-q58UM-bsgyq9ScHa4hfNOWGelM",
    "location": {
        "digest": "XHzVIQ1JEps1jDJIMGF_kn7YIcQdn5bW",
        "id": "ga4gh:SL.XHzVIQ1JEps1jDJIMGF_kn7YIcQdn5bW",
        "end": 20003010,
        "start": 20003010,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "G", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1264314/?new_evidence=true
duplication_inputs = {
    "hgvs": "NC_000013.11:g.19993838_19993839dup",
    "spdi": "NC_000013.11:19993837:GT:GTGT",
    "gnomad": "13-19993838-GT-GTGT",
}

duplication_output = {
    "digest": "2Ju2sXpBNposOefRvorsmfqkAvt9tRHD",
    "id": "ga4gh:VA.2Ju2sXpBNposOefRvorsmfqkAvt9tRHD",
    "location": {
        "digest": "yCVGYQzbSLQe-GeAaHbW0dOiEGzHF3Yj",
        "id": "ga4gh:SL.yCVGYQzbSLQe-GeAaHbW0dOiEGzHF3Yj",
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {
        "length": 4,
        "repeatSubunitLength": 2,
        "sequence": "GTGT",
        "type": "ReferenceLengthExpression",
    },
    "type": "Allele",
}


def test_normalized_results(caching_translator):
    tlr = caching_translator
    tlr.normalize = True

    assert (
        tlr._from_gnomad(snv_inputs["gnomad"]).model_dump(exclude_none=True)
        == snv_output
    )
    assert (
        tlr._from_gnomad(deletion_inputs["gnomad"]).model_dump(exclude_none=True)
        == gnomad_deletion_output
    )
    assert (
        tlr._from_gnomad(insertion_inputs["gnomad"]).model_dump(exclude_none=True)
        == gnomad_insertion_output
    )
    assert (
        tlr._from_gnomad(duplication_inputs["gnomad"]).model_dump(exclude_none=True)
        == duplication_output
    )


def test_results(caching_translator):
    """Ensure we can get the same results from gnomad as vrs-python, ie
    that the id and digest are computed recursively for the Allele object"""
    tlr = caching_translator
    assert tlr is not None

    inputs_dict = {
        "snv": (snv_inputs, snv_output),
        "deletion": (deletion_inputs, gnomad_deletion_output),
        "insertion": (insertion_inputs, gnomad_insertion_output),
        "duplication": (duplication_inputs, duplication_output),
    }

    for variant_type, (input, expected_allele) in inputs_dict.items():
        gnomad_expr = input["gnomad"]
        # allele object validation
        allele_id = tlr.translate_from(fmt="gnomad", var=gnomad_expr)

        assert expected_allele["id"] == allele_id, (
            f"{allele_id} does not match for {gnomad_expr}: "
            + f"\nexpected: {expected_allele["id"]}"
        )


def test_cache(caching_translator):
    """Ensure that results from getting allele IDs are faster the second time."""
    tlr = caching_translator
    assert tlr is not None
    all_inputs = [snv_inputs, deletion_inputs, insertion_inputs, duplication_inputs]

    start_time = time.time()
    for inputs in all_inputs:
        tlr.translate_from(fmt="gnomad", var=inputs["gnomad"])
    noncached_time = time.time() - start_time

    start_time = time.time()
    for inputs in all_inputs:
        tlr.translate_from(fmt="gnomad", var=inputs["gnomad"])
    cache_time = time.time() - start_time

    assert cache_time < (
        noncached_time / 4
    ), f"Cache should make things significantly faster first {noncached_time} second {cache_time}."


@pytest.fixture()
def num_threads():
    """Return the number of threads to use for testing."""
    # TODO - anything more than 2 produces errors, need to investigate
    # e.g.  {'error': 'bad parameter or other API misuse', 'item': {'fmt': 'gnomad', 'var': '13-20003010-A-AG'}}
    return 20


def test_threading(translator, num_threads):
    """Ensure we can feed the threaded_translate_from method a generator and get results back."""
    tlr = translator
    assert tlr is not None

    parameters = [
        {"fmt": "gnomad", "var": snv_inputs["gnomad"]},
        {"fmt": "gnomad", "var": deletion_inputs["gnomad"]},
        {"fmt": "gnomad", "var": insertion_inputs["gnomad"]},
        {"fmt": "gnomad", "var": duplication_inputs["gnomad"]},
    ]

    def repeat_sequence(sequence, times):
        """
        Generator that repeats a given sequence for a specified number of times.

        Args:
        - sequence: The sequence to be repeated.
        - times: The number of times to repeat the sequence.

        Yields:
        - Elements from the sequence in a repeated pattern.
        """
        line_number = 1
        for _ in range(times):
            for _sequence in sequence:
                yield VCFItem(
                    **_sequence,
                    file_name="test",
                    line_number=line_number,
                    identifier=None,
                )
                line_number += 1

    c = 0  # count of results
    _times = 2
    for result in tlr.translate_from(
        repeat_sequence(parameters, times=_times), num_threads=num_threads
    ):
        c += 1
        validate_threaded_result(result, validate_passthrough=False)


@pytest.fixture
def gnomad_csv() -> pathlib.Path:
    """Return a path to a gnomad csv file."""
    _ = pathlib.Path(
        "tests/fixtures/gnomAD_v4.0.0_ENSG00000012048_2024_03_04_18_33_26.csv"
    )
    assert _.exists()
    return _


def gnomad_ids(path, limit=None) -> Generator[tuple, None, None]:
    """Open the csv file and yield the first column 'gnomAD ID', skip the header"""
    c = 0
    with open(path, "r") as f:
        skip = True
        for line in f:
            if skip:
                skip = False
                continue
            _ = line.split(",")
            gnomad_id = _[0]
            # use allele number as the identifier
            allele_number = _[18]
            yield VCFItem("gnomad", gnomad_id, path, c, allele_number)
            c += 1
            if limit and c > limit:
                break


def test_gnomad(translator, gnomad_csv, num_threads):
    """We can process a set of gnomad variants."""
    tlr = translator
    assert tlr is not None
    c = 0  # count of results
    for result_dict in tlr.translate_from(
        generator=gnomad_ids(gnomad_csv), num_threads=num_threads
    ):
        c += 1
        validate_threaded_result(result_dict, validate_passthrough=False)


def test_gnomad_inline(translator, gnomad_csv, num_threads=1):
    """We can process a set of gnomad variants."""
    tlr = translator
    assert tlr is not None
    c = 0  # count of results
    for result_dict in tlr.translate_from(
        generator=gnomad_ids(gnomad_csv), num_threads=num_threads
    ):
        c += 1
        validate_threaded_result(result_dict, validate_passthrough=False)
