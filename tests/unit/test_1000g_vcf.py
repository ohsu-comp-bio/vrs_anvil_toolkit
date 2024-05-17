from tests.unit import validate_threaded_result
from vrs_anvil import params_from_vcf


def test_1000g_vcf(translator, thousand_genome_vcf, num_threads):
    """Read and process AnVIL sourced vcf."""
    tlr = translator
    assert tlr is not None
    c = 0  # count of results
    for result_dict in tlr.translate_from(
        generator=params_from_vcf(thousand_genome_vcf), num_threads=num_threads
    ):
        c += 1
        print(result_dict)
        validate_threaded_result(result_dict, validate_passthrough=True)

    validate_threaded_result(result_dict)
