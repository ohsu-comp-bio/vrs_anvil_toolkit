import pathlib
import time
from typing import Generator

import pytest

# see https://github.com/ga4gh/vrs-python/blob/main/tests/extras/test_allele_translator.py#L17

snv_inputs = {
    "hgvs": "NC_000019.10:g.44908822C>T",
    "beacon": "19 : 44908822 C > T",
    "spdi": "NC_000019.10:44908821:1:T",
    "gnomad": "19-44908822-C-T"
}

snv_output = {
    "location": {
        "end": 44908822,
        "start": 44908821,
        "sequenceReference": {
            "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "T",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1373966/?new_evidence=true
deletion_inputs = {
    "hgvs": "NC_000013.11:g.20003097del",
    "spdi": ["NC_000013.11:20003096:C:", "NC_000013.11:20003096:1:"],
    "gnomad": "13-20003096-AC-A"
}

deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

gnomad_deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003095,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "A",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

deletion_output_normalized = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "length": 0,
        "repeatSubunitLength": 1,
        "sequence": "",
        "type": "ReferenceLengthExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1687427/?new_evidence=true
insertion_inputs = {
    "hgvs": "NC_000013.11:g.20003010_20003011insG",
    "spdi": ["NC_000013.11:20003010::G", "NC_000013.11:20003010:0:G"],
    "gnomad": "13-20003010-A-AG"
}

insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003010,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "G",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

gnomad_insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003009,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "AG",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1264314/?new_evidence=true
duplication_inputs = {
    "hgvs": "NC_000013.11:g.19993838_19993839dup",
    "spdi": "NC_000013.11:19993837:GT:GTGT",
    "gnomad": "13-19993838-GT-GTGT"
}

duplication_output = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "GTGT",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

duplication_output_normalized = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "length": 4,
        "repeatSubunitLength": 2,
        "sequence": "GTGT",
        "type": "ReferenceLengthExpression"
    },
    "type": "Allele"
}


def _remove_ids(_dict):
    """Remove the "id" key from a dictionary TODO - why is this necessary?"""
    if "id" in _dict:
        del _dict["id"]
    if "location" in _dict and "id" in _dict["location"]:
        del _dict["location"]["id"]
    return _dict


def test_cache(my_translator):
    """Ensure we can get the same results from gnomad as vrs-python and that results are faster the second time."""
    tlr = my_translator
    assert tlr is not None
    tlr.normalize = False

    start_time = time.time()

    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=snv_inputs["gnomad"]).model_dump(exclude_none=True)) == snv_output
    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=deletion_inputs["gnomad"]).model_dump(exclude_none=True)) == gnomad_deletion_output
    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=insertion_inputs["gnomad"]).model_dump(exclude_none=True)) == gnomad_insertion_output
    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=duplication_inputs["gnomad"]).model_dump(exclude_none=True)) == duplication_output

    end_time = time.time()

    elapsed_time = end_time - start_time
    # print(f"Elapsed time: {elapsed_time} seconds")

    start_time = time.time()

    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=snv_inputs["gnomad"]).model_dump(exclude_none=True)) == snv_output
    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=deletion_inputs["gnomad"]).model_dump(exclude_none=True)) == gnomad_deletion_output
    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=insertion_inputs["gnomad"]).model_dump(exclude_none=True)) == gnomad_insertion_output
    assert _remove_ids(tlr.translate_from(fmt="gnomad", var=duplication_inputs["gnomad"]).model_dump(exclude_none=True)) == duplication_output

    end_time = time.time()

    elapsed_time2 = end_time - start_time
    # print(f"Elapsed time: {elapsed_time2} seconds")

    assert elapsed_time2 < (elapsed_time / 100), f"Cache should make things significantly faster first {elapsed_time} second {elapsed_time2}."


@pytest.fixture()
def num_threads():
    """Return the number of threads to use for testing."""
    # TODO - anything more than 2 produces errors, need to investigate
    # e.g.  {'error': 'bad parameter or other API misuse', 'item': {'fmt': 'gnomad', 'var': '13-20003010-A-AG'}}
    return 20


def test_threading(threaded_translator, num_threads):
    """Ensure we can feed the threaded_translate_from method a generator and get results back."""
    tlr = threaded_translator
    assert tlr is not None
    tlr.normalize = False

    parameters = [
        {"fmt": "gnomad", "var": snv_inputs["gnomad"]},
        {"fmt": "gnomad", "var": deletion_inputs["gnomad"]},
        {"fmt": "gnomad", "var": insertion_inputs["gnomad"]},
        {"fmt": "gnomad", "var": duplication_inputs["gnomad"]}
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
        for _ in range(times):
            yield from sequence

    c = 0  # count of results
    _times = 2
    errors = []
    for result in tlr.threaded_translate_from(repeat_sequence(parameters, times=_times), num_threads=num_threads):
        c += 1
        assert result is not None
        if isinstance(result, dict):
            _ = result
        else:
            _ = result.model_dump(exclude_none=True)
        assert _ is not None
        if 'error' in _:
            errors.append(_)
            continue
        for k in ['location', 'state', 'type']:
            assert k in _
    assert c == (_times * len(parameters)), f"Expected {_times * len(parameters)} results, got {c}."
    assert len(errors) == 0, f"num_threads {num_threads} {c} items {len(errors)} errors {errors}."


@pytest.fixture
def gnomad_csv() -> pathlib.Path:
    """Return a path to a gnomad csv file."""
    _ = pathlib.Path("tests/data/gnomAD_v4.0.0_ENSG00000012048_2024_03_04_18_33_26.csv")
    assert _.exists()
    return _


def gnomad_ids(path, limit=None) -> Generator[dict, None, None]:
    """Open the csv file and yield the first column 'gnomAD ID', skip the header"""
    c = 0
    with open(path, "r") as f:
        skip = True
        for line in f:
            if skip:
                skip = False
                continue
            gnomad_id = line.split(",")[0]
            yield {"fmt": "gnomad", "var": gnomad_id}
            c += 1
            if limit and c > limit:
                break


def test_gnomad(threaded_translator, gnomad_csv, num_threads):
    """We can process a set of gnomad variants."""
    tlr = threaded_translator
    assert tlr is not None
    c = 0  # count of results
    start_time = time.time()
    errors = []
    for result in tlr.threaded_translate_from(generator=gnomad_ids(gnomad_csv), num_threads=num_threads):
        c += 1
        assert result is not None
        if isinstance(result, dict):
            _ = result  # handle dummy results
        else:
            _ = result.model_dump(exclude_none=True)
        assert _ is not None
        if 'error' in _:
            errors.append(_)
            continue
        for k in ['location', 'state', 'type']:
            assert k in _
    end_time = time.time()

    elapsed_time = end_time - start_time

    print(errors)
    assert len(errors) == 0, f"num_threads {num_threads} elapsed time: {elapsed_time} seconds {c} items {len(errors)} errors {errors}."
