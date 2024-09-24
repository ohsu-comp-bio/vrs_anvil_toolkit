import os
from typing import Generator

import pytest
from pysam import VariantFile, VariantRecord
from datetime import datetime


@pytest.fixture
def chromosome():
    """Return chromosome in the VCF."""
    return "chrY"


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
    assert "GT" in record.format, "Genotype (GT) is required"

    for participant, values in record.samples.items():
        assert "GT" in values, "Genotype (GT) is required"
        # see https://samtools.github.io/hts-specs/VCFv4.1.pdf

        # TODO - this test is a bit naive, should we be more robust. consider AD, GQ, RGQ?
        if any(values["GT"]):
            yield participant


def test_remote_vcf(vcf_path, chromosome, start, stop, expected_record_count):
    """Read a remote vcf file, query a range of alleles, check that at least 1 participant exists for each allele."""
    assert "GCS_OAUTH_TOKEN" in os.environ, (
        "GCS_OAUTH_TOKEN is required "
        "export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token) "
        "see https://github.com/pysam-developers/pysam/issues/592#issuecomment-353693674 https://support.terra.bio/hc/en-us/articles/360042647612-May-4-2020"
    )
    try:
        vcf = VariantFile(vcf_path)  # auto-detect input format
        # fetch returns pysam.libcbcf.VariantRecord
        records = [_ for _ in vcf.fetch(chromosome, start, stop)]
        assert len(records) == expected_record_count

        for variant_record in records:
            my_participants = [_ for _ in participants(variant_record)]
            assert len(my_participants) < len(
                variant_record.samples
            ), "Not all participants have this allele"
            assert len(my_participants) > 0, "No participants have this allele"
    except ValueError as e:
        print("ValueError: has GCS_OAUTH_TOKEN expired?", e)
        raise e


def test_caf_generation_given_row(vcf_path):
    vcf = VariantFile(vcf_path)

    for record in vcf.fetch():
        # FIXME: pulling out vrs ids has to mirror the params of the VRS annotator...
        # eg: if compute_for_ref = True, then
        vrs_ids = record.info["VRS_Allele_IDs"]
        print("VRS IDs:", vrs_ids)

        focus_allele_count = 0
        locus_allele_count = 0
        num_homozygotes = 0
        num_hemizygotes = 0

        assert len(vcf.header.samples) == len(
            record.samples
        ), "row has different number of samples as the record"

        # print([a.allele_indices for a in record.samples.values()])
        for patient_id, genotype in record.samples.items():
            alleles = genotype.allele_indices

            # less than total number of samples, as only samples with a variant read are counted

            # increment counts only if variant read/measured in cohort
            if alleles == (None, None):
                continue

            num_alleles = sum(alleles)

            # record zygosity
            if num_alleles == 1:
                num_hemizygotes += 1
            elif num_alleles == 2:
                num_homozygotes

            # increment counts
            focus_allele_count += num_alleles
            locus_allele_count += 2

        # populate caf dict
        # FIXME: how do we figure out actual VRS ID if multiple?
        allele_id = vrs_ids[-1]
        allele_frequency = focus_allele_count * 1.0 / locus_allele_count

        caf_dict = {
            "type": "CohortAlleleFrequency",
            "label": f"Overall Cohort Allele Frequency for {allele_id}",
            "derivedFrom": {
                "id": "GREGOR_COMBINED_CONSORTIUM_U07",
                "type": "DataSet",
                "label": "GREGoR Combined Consortium U07",
                "version": f"Created {datetime.now()}",
            },
            "focusAllele": allele_id,
            "focusAlleleCount": focus_allele_count,
            "locusAlleleCount": locus_allele_count,
            "alleleFrequency": allele_frequency,
            "cohort": {
                "id": "GREGOR_COMBINED_CONSORTIUM_U07",
                "label": "GREGoR Combined Consortium U07",
            },
            "ancillaryResults": {
                "homozygotes": num_homozygotes,
                "hemizygotes": num_hemizygotes,
            },
        }
        import json

        print(json.dumps(caf_dict))
        break

    assert caf_dict["focusAlleleCount"] == 10, "expected focusAlleleCount of 10"
    assert caf_dict["locusAlleleCount"] == 1168, "expected locusAlleleCount of 1168"
