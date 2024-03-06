

def generate_gnomad_id(vcf_line):
    # Assuming a standard VCF format with tab-separated fields
    fields = vcf_line.strip().split('\t')

    # Extract relevant information (you may need to adjust these indices based on your VCF format)
    chromosome = fields[0]
    position = fields[1]
    reference_allele = fields[3]
    alternate_allele = fields[4]

    # Create a gnomAD-like ID
    gnomad_id = f"{chromosome}-{position}-{reference_allele}-{alternate_allele}"

    return gnomad_id


def test_generate_gnomad_id(my_translator):
    input_vcf = f"tests/data/test_vcf_input.vcf"
    errors = []
    results = []
    with open(input_vcf, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            gnomad_id = generate_gnomad_id(line)
            try:
                result = my_translator.translate_from(fmt="gnomad", var=gnomad_id)
                results.append((gnomad_id, result))
            except Exception as e:
                errors.append((gnomad_id, e))
    print(errors)
    assert len(errors) == 0, f"Errors: {len(errors)} Successes: {len(results)}"
