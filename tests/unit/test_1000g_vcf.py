import time

from tests.unit import validate_threaded_result
from vrs_anvil import params_from_vcf


def test_1000g_vcf(threaded_translator, thousand_genome_vcf, num_threads):
    """Read and process AnVIL sourced vcf."""
    tlr = threaded_translator
    assert tlr is not None
    c = 0  # count of results
    start_time = time.time()
    errors = []
    for result_dict in tlr.threaded_translate_from(
        generator=params_from_vcf(thousand_genome_vcf), num_threads=num_threads
    ):
        c += 1
        validate_threaded_result(result_dict, errors, validate_passthrough=True)
    end_time = time.time()

    elapsed_time = end_time - start_time

    print(errors)
    assert (
        len(errors) == 0
    ), f"num_threads {num_threads}, elapsed time: {elapsed_time} seconds, {c} items, {len(errors)} errors {errors}."
