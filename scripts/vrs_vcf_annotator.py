from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator
from vrs_anvil import seqrepo_dir

SEQREPO_DIR = seqrepo_dir()
input_vcf = "tests/fixtures/1kGP.chr1.1000.vcf"
input_vcf_path_stem = input_vcf.replace(".vcf", "")
output_pkl = f"{input_vcf_path_stem}-test-vrs-objects.pkl"

output_vcf = None
require_validation = True
rle_seq_limit = 50

vcf_annotator = VCFAnnotator(seqrepo_root_dir=SEQREPO_DIR)
vcf_annotator.tlr.rle_seq_limit = rle_seq_limit
vcf_annotator.annotate(
    vcf_in=input_vcf,
    vcf_out=output_vcf,
    vrs_pickle_out=output_pkl,
    require_validation=require_validation,
)

print(f"written to {output_pkl}")
