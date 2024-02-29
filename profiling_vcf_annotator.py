from utils import annotate_vcf, seqrepo_dir

SEQREPO_DIR = seqrepo_dir()
input_vcf = "vcf/NA12878.chr1.100.000.vcf"
stem = input_vcf.replace('.vcf', '')
# output_vcf = f"{stem}output.vcf.gz"
output_pkl = f"{stem}-vrs-objects.pkl"
output_vcf = None
# output_pkl = None

annotate_vcf(input_vcf, output_vcf, output_pkl, SEQREPO_DIR)