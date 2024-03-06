import logging

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)


def generate_gnomad_ids(vcf_line, compute_for_ref: bool = True) -> list[str]:
    # Assuming a standard VCF format with tab-separated fields
    fields = vcf_line.strip().split('\t')
    gnomad_ids = []
    # Extract relevant information (you may need to adjust these indices based on your VCF format)
    chromosome = fields[0]
    position = fields[1]
    reference_allele = fields[3]
    alternate_allele = fields[4]

    gnomad_loc = f"{chromosome}-{position}"
    if compute_for_ref:
        gnomad_ids.append(f"{gnomad_loc}-{reference_allele}-{reference_allele}")
    for alt in alternate_allele.split(","):
        alt = alt.strip()
        if '*' in alt:
            _logger.debug("Star allele found: %s", alt)
            continue
        gnomad_ids.append(f"{gnomad_loc}-{reference_allele}-{alt}")

    return gnomad_ids


def test_generate_gnomad_ids(my_translator):
    input_vcf = f"tests/data/test_vcf_input.vcf"
    errors = []
    results = []
    with open(input_vcf, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            for gnomad_id in generate_gnomad_ids(line):
                try:
                    result = my_translator.translate_from(fmt="gnomad", var=gnomad_id)
                    results.append((gnomad_id, result))
                except Exception as e:
                    errors.append((gnomad_id, e))
    print(errors)
    assert len(results) >= 12, f"Errors: {len(errors)} Successes: {len(results)}"
    # TODO confirm these are expected errors ? [('chr19-54220999-A-A', ValidationError('Expected reference sequence A on GRCh38:chr19 at positions (54220998, 54220999) but found T')), ('chr19-54220999-A-T', ValidationError('Expected reference sequence A on GRCh38:chr19 at positions (54220998, 54220999) but found T')), ('chr19-54221654-A-A', ValidationError('Expected reference sequence A on GRCh38:chr19 at positions (54221653, 54221654) but found T')), ('chr19-54221654-A-T', ValidationError('Expected reference sequence A on GRCh38:chr19 at positions (54221653, 54221654) but found T')), ('chr19-54221654-A-P', ValueError('Unable to parse data as gnomad variation'))]
