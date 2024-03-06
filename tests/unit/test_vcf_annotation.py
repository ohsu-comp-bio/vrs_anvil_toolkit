import os
import pytest
import vcf

from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator

VRS_ERROR_FIELD = VCFAnnotator.VRS_ERROR_FIELD
VRS_ALLELE_IDS_FIELD = VCFAnnotator.VRS_ALLELE_IDS_FIELD
TEST_DATA_DIR = "tests/data"

@pytest.fixture
def seqrepo_dir():
   with open(".env") as f:
      for line in f:
         # Ignore comments and empty lines
         if line.strip() and not line.strip().startswith('#'):
            key, value = line.strip().split('=', 1)
            if key == "SEQREPO_ROOT":
               return value + "/latest"

@pytest.fixture(scope='function')
def delete_output(request):
   # cleanup output file after test execution
   def cleanup():
      os.remove(request.node.output_vcf)

   request.addfinalizer(cleanup)


@pytest.mark.skip(reason="Issue with test data see https://github.com/ohsu-comp-bio/vrs-python-testing/issues/16")
def test_small_vcf_annotation(seqrepo_dir, delete_output, request):
   
   # create annotated vcf test file 
   input_vcf = f"{TEST_DATA_DIR}/test_vcf_input.vcf"
   expected_vcf = f"{TEST_DATA_DIR}/expected_vcf_output.vcf"
   output_vcf = f"{TEST_DATA_DIR}/test_vcf_output.vcf"
   request.node.output_vcf = output_vcf

   vcf_annotator = VCFAnnotator(seqrepo_root_dir=seqrepo_dir) 
   vcf_annotator.annotate(input_vcf, output_vcf)

   # assert ids are the same line by line
   test_vcf_reader = vcf.Reader(open(output_vcf, 'r'))
   expected_vcf_reader = vcf.Reader(open(expected_vcf, 'r'))
   
   for i, expected_record in enumerate(expected_vcf_reader):
      try:
         test_record = next(test_vcf_reader)
      except:
         assert False, "missing variants, not enough rows to compare to"
         
      variant_str = f"{i+1}th variant with record {test_record}"

      # check successful annotation
      # nice to have: output entire error, as pyvcf delimits with space
      assert VRS_ERROR_FIELD not in test_record.INFO, \
         f"VRS ERROR: unable to annotate {variant_str}, please check it"
      
      # check each ID matches
      test_ids = test_record.INFO[VRS_ALLELE_IDS_FIELD]
      expected_ids = expected_record.INFO[VRS_ALLELE_IDS_FIELD]

      if isinstance(expected_ids, list):
         assert isinstance(test_ids, list) and len(expected_ids) == len(test_ids), \
            f"not the same number of IDs for {variant_str}"

         for expected_id, test_id in zip(expected_ids, test_ids):
            assert expected_id == test_id, \
               f"ids for {variant_str} do not match"