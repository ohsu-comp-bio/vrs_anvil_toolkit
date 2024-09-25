import json
import os
import pandas as pd
import pytest

from datetime import datetime
from typing import Generator
from pysam import VariantFile, VariantRecord


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


def test_caf_generation_from_vcf_row(vcf_path):
    vcf = VariantFile(vcf_path)
    ref_computed = True
    caf_dicts = []

    for i, record in enumerate(vcf.fetch()):
        vrs_ids = record.info["VRS_Allele_IDs"]
        print("VRS IDs:", vrs_ids)
        print(len(vrs_ids) - 1)

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
                num_homozygotes += 1

            # increment counts
            focus_allele_count += num_alleles
            locus_allele_count += 2

        allele_frequency = focus_allele_count * 1.0 / locus_allele_count

        # populate caf dict(s) with each allele, multiple if multiple alts per variant
        if ref_computed:
            allele_ids = vrs_ids[1:]
        else:
            allele_ids = vrs_ids

        for allele_id in allele_ids:
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

            print(json.dumps(caf_dict))
            caf_dicts.append(caf_dict)

        # check first caf
        if i == 0:
            assert caf_dict["focusAlleleCount"] == 10, "expected focusAlleleCount of 10"
            assert (
                caf_dict["locusAlleleCount"] == 1168
            ), "expected locusAlleleCount of 1168"

        # do only the first 10
        if i == 9:
            break

    expected_num_cafs = 13
    assert (
        len(caf_dicts) == expected_num_cafs
    ), f"expected {expected_num_cafs} caf_dicts, got {len(caf_dicts)}"


def test_get_phenotypes_from_vcf_row(vcf_path):
    vcf = VariantFile(vcf_path)

    for i, record in enumerate(vcf.fetch()):
        vrs_ids = record.info["VRS_Allele_IDs"]
        phenotypes_set = set()
        patients = [
            patient
            for patient, genotype in record.samples.items()
            if 1 in genotype.allele_indices
        ]
        print(f"{record.chrom}-{record.pos}-{record.ref}-{record.alts}")

        print("VRS IDs:", vrs_ids)

        for patient_id in patients:
            phenotypes = get_patient_phenotypes(
                patient_id, "/Users/wongq/data/gregor/phenotypes.tsv"
            )
            phenotypes_set.update(phenotypes)

        print(phenotypes_set)
        print("len(phenotypes_set):", len(phenotypes_set))

        if i == 3:
            expected_num_phenotypes = 225
            assert (
                len(phenotypes_set) == expected_num_phenotypes
            ), f"Expected {expected_num_phenotypes} phenotypes, got {len(phenotypes_set)}"

        if i == 9:
            break


# FIXME: general methods to also move out of tests
def get_patient_phenotypes(
    patient_id: str, phenotypes_table: str = None, cached_dict: str = None
) -> list:
    """get list of phenotypes associated with a patient"""

    # load from cache if exists
    if cached_dict is not None and os.path.exists(cached_dict):
        with open(cached_dict, "r") as file:
            patient_to_phenotypes = json.load(file)
            return patient_to_phenotypes[patient_id]

    # if unspecified, parse table using Terra table
    if phenotypes_table is None:
        assert False, "pulling from Terra phenotypes table not implemented yet"

    # table path specified, parse using that table
    with open(phenotypes_table, "r") as file:
        if phenotypes_table.endswith(".csv"):
            pheno_df = pd.read_csv(file)
        elif phenotypes_table.endswith(".tsv"):
            pheno_df = pd.read_csv(file, sep="\t")
        else:
            raise Exception(
                "Only csv and tsv file types implemented for phenotype table"
            )

        phenotypes = pheno_df[pheno_df["participant_id"] == patient_id][
            "term_id"
        ].unique()

        return list(phenotypes)
