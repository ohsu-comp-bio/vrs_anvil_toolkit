import filecmp
import os
import pytest

from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator

TEST_DATA_DIR = "tests/data"

def get_latest_seqrepo_dir():
    with open(".env") as f:
        for line in f:
            # Ignore comments and empty lines
            if line.strip() and not line.strip().startswith('#'):
                key, value = line.strip().split('=', 1)
                if key == "SEQREPO_LATEST_DIR":
                    return value

def test_small_vcf_annotation():
    input_vcf = f"{TEST_DATA_DIR}/test_vcf_input.vcf"
    output_vcf = f"{TEST_DATA_DIR}/test_vcf_output.vcf.gz"
    output_vrs_pkl = f"{TEST_DATA_DIR}/test_vrs_objects.pkl"
    expected_vrs_pkl = f"{TEST_DATA_DIR}/expected_vrs_objects.pkl"

    vcf_annotator = VCFAnnotator(seqrepo_root_dir=get_latest_seqrepo_dir())
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl)
    
    try:
        assert filecmp.cmp(output_vrs_pkl, expected_vrs_pkl), "files not identical"
    finally:
        os.remove(output_vcf)
        os.remove(output_vrs_pkl)