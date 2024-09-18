import os
from typing import Generator

import pytest
from pysam import VariantFile, VariantRecord


@pytest.fixture
def url():
    """Return a remote VCF URL."""
    return os.environ['MY_OBJECT']


@pytest.fixture
def chromosome():
    """Return chromosome in the VCF."""
    return 'chrY'


@pytest.fixture
def start():
    """Return the start range to query."""
    return 2786855


@pytest.fixture
def stop():
    """Return the end range to query."""
    return 2787682


@pytest.fixture
def expected_record_count():
    """Return the expected record count."""
    return 2


def participants(record: VariantRecord) -> Generator[str, None, None]:
    """Return the participants that `have` this allele."""
    assert 'GT' in record.format, "Genotype (GT) is required"

    for participant, values in record.samples.items():
        assert 'GT' in values, "Genotype (GT) is required"
        # see https://samtools.github.io/hts-specs/VCFv4.1.pdf

        # In Variant Call Format (VCF) files, GT stands for genotype, which is encoded as allele values separated by a slash (/) or vertical pipe (|):
        # 0: The reference base
        # 1: The first entry in the ALT column
        # 2: The second allele listed in ALT
        # Forward slash (/): Indicates that no phasing information is available
        # Vertical pipe (|): Indicates that the genotype is phased

        # TODO - this test is a bit naive, should we be more robust. consider AD, GQ, RGQ?
        if any(values['GT']):
            yield participant


        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">

        # see https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format#:~:text=Allele%20depth%20(AD)%20and%20depth,each%20of%20the%20reported%20alleles.
        # In a Variant Call Format (VCF) file, AD stands for Allele Depth, which is the number of reads that support each allele. The AD field is an array that includes the reference allele as the first entry, and the remaining entries are for each alternate allele at that locus

        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">

        # see https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format#:~:text=PL%20is%20calculated.-,GQ,how%20they%20should%20be%20used.
        # In Variant Call Format (VCF), GQ stands for Genotype Quality and is a measure of how confident a genotype assignment is. GQ is calculated by taking the difference between the PL of the second most likely genotype and the PL of the most likely genotype. The PLs are normalized so that the most likely PL is always 0, so the GQ is usually equal to the second smallest PL. However, the GQ is capped at 99, so if the second most likely PL is greater than 99, the GQ will be 99.

        ##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">

        # see https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized
        # Reference Genotype Quality (RGQ) -- The phred-scaled confidence that the reference genotypes are correct.  A higher score indicates a higher confidence.  For more information on RGQ, please see the GQ documentation, but note that RGQ applies to the reference, not the variant.  For more information on interpreting phred-scaled values, please see Phred-scaled quality scores.


def test_remote_vcf(url, chromosome, start, stop, expected_record_count):
    """Read a remote vcf file, query a range of alleles, check that at least 1 participant exists for each allele."""
    assert 'GCS_OAUTH_TOKEN' in os.environ, "GCS_OAUTH_TOKEN is required\nexport GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)\nsee https://github.com/pysam-developers/pysam/issues/592#issuecomment-353693674 https://support.terra.bio/hc/en-us/articles/360042647612-May-4-2020"
    try:
        vcf = VariantFile(url)  # auto-detect input format
        # fetch returns pysam.libcbcf.VariantRecord
        records = [_ for _ in vcf.fetch(chromosome, start, stop)]
        assert len(records) == expected_record_count
        for variant_record in records:
            my_participants = [_ for _ in participants(variant_record)]
            assert len(my_participants) < len(variant_record.samples), "Not all participants have this allele"
            assert len(my_participants) > 0, "No participants have this allele"
    except ValueError as e:
        print("ValueError: has GCS_OAUTH_TOKEN expired?", e)
        raise e
