import json
import os
import pytest

from typing import Generator
from pysam import VariantFile, VariantRecord

from vrs_anvil.phenotype import get_cohort_allele_frequency, get_patient_phenotype_index

############
# FIXTURES #
############


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


@pytest.fixture
def remote_vcf_path():
    return "gs://fc-secure-0392c60d-b346-4e60-b5f1-602a9bbfb0e1/gregor_joint_callset_022824/chr_vcfs/gregor_consortium_dataset_u03_v2.chrY.vcf.gz"


@pytest.fixture()
def vrs_id(vcf_path):
    for i, record in enumerate(VariantFile(vcf_path)):
        if i == 3:
            return record.info["VRS_Allele_IDs"][1]  # 1 since 0 is a ref


@pytest.fixture()
def phenotype_table():
    assert (
        "PHENOTYPE_TABLE" in os.environ
    ), "no PHENOTYPE_TABLE bash variable defining the path of the phenotype table, make sure to export"

    return os.environ["PHENOTYPE_TABLE"]


def participants(record: VariantRecord) -> Generator[str, None, None]:
    """Return the participants that `have` this allele."""
    assert "GT" in record.format, "Genotype (GT) is required"

    for participant, values in record.samples.items():
        assert "GT" in values, "Genotype (GT) is required"
        # see https://samtools.github.io/hts-specs/VCFv4.1.pdf

        # TODO - this test is a bit naive, should we be more robust. consider AD, GQ, RGQ?
        if any(values["GT"]):
            yield participant


def test_remote_vcf(remote_vcf_path, chromosome, start, stop, expected_record_count):
    """Read a remote vcf file, query a range of alleles, check that at least 1 participant exists for each allele."""
    assert "GCS_OAUTH_TOKEN" in os.environ, (
        "GCS_OAUTH_TOKEN required: "
        "export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)"
        "see https://github.com/pysam-developers/pysam/issues/592#issuecomment-353693674 https://support.terra.bio/hc/en-us/articles/360042647612-May-4-2020"
    )
    try:
        vcf = VariantFile(remote_vcf_path)  # auto-detect input format
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


def test_get_phenotypes_from_vcf_row(vcf_path, phenotype_table):
    """calculate list of unique phenotypes for the first 4 rows, checking total unique phenotypes the last row"""

    vcf = VariantFile(vcf_path)
    assert "VRS_Allele_IDs" in vcf.header.info, (
        "no VRS_Allele_IDs key in INFO found, "
        "please ensure that this is an VRS annotated VCF"
    )

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

        phenotype_index = get_patient_phenotype_index(phenotype_table)

        for patient_id in patients:
            if patient_id in phenotype_index:
                phenotypes_set.update(phenotype_index[patient_id])
        print("len(phenotypes_set):", len(phenotypes_set))

        if i == 3:
            expected_num_phenotypes = 225
            assert (
                len(phenotypes_set) == expected_num_phenotypes
            ), f"Expected {expected_num_phenotypes} phenotypes, got {len(phenotypes_set)}"

            break


def test_allele_counts_first_10_rows(phenotype_table):
    """test that the calculated allele counts with no phenotype specified matches the actual counts stored in first 10 rows. Works for diploid (non-sex) variants"""

    vcf_path = "/Users/wongq/data/gregor/ann.100k.chr3.vcf.gz"
    vcf = VariantFile(vcf_path)
    has_ref = "REF" in vcf.header.info["VRS_Allele_IDs"].description

    for i, record in enumerate(vcf.fetch()):
        print("~~~~~~~~~ row ", i, "~~~~~~~~~")

        # use only alt VRS IDs, not REFs
        vrs_allele_ids = record.info["VRS_Allele_IDs"]
        print("vrs_allele_ids:", vrs_allele_ids)
        if has_ref:
            vrs_allele_ids = vrs_allele_ids[1:]

        # for each alt ID, ensure stored allele counts match calculated allele counts
        for alt_index, allele_id in enumerate(vrs_allele_ids):
            print("alt id:", allele_id)
            if not allele_id:
                continue

            caf = get_cohort_allele_frequency(
                allele_id, vcf_path, phenotype_table=phenotype_table
            )

            print("alt_index", alt_index)
            print("AC:", record.info["AC"][alt_index], caf["focusAlleleCount"])
            print("AN:", record.info["AN"], caf["locusAlleleCount"])
            ac = record.info["AC"][alt_index]
            an = record.info["AN"]

            assert (
                caf["focusAlleleCount"] == ac
            ), f"row {i} alt {alt_index} has different focus allele counts, expected {ac}, got {caf['focusAlleleCount']}"
            assert (
                caf["locusAlleleCount"] == an
            ), f"row {i} alt {alt_index} has different locus allele counts, expected {an}, got {caf['locusAlleleCount']}"

        if i == 10:
            break


def test_default_caf_gen(vrs_id, vcf_path, phenotype_table):
    """test caf generation with default parameters and no phenotype specified"""
    caf = get_cohort_allele_frequency(vrs_id, vcf_path, phenotype_table=phenotype_table)
    print(json.dumps(caf))

    # sanity checks
    assert (
        caf["type"] == "CohortAlleleFrequency"
    ), f"object of type CohortAlleleFrequency not returned, {caf['type']} instead"
    assert (
        caf["focusAlleleCount"] <= caf["locusAlleleCount"]
    ), f"Focus allele count ({caf['focusAlleleCount']}) is larger than locus allele count ({caf['locusAlleleCount']})"

    print("focusAlleleCount:", caf["focusAlleleCount"])
    print("locusAlleleCount:", caf["locusAlleleCount"])

    # check allele frequency
    expected_allele_freq = 0.9465
    actual_allele_freq = round(caf["alleleFrequency"], 4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"

    # ensure fields exist
    expected_fields = [
        "focusAlleleCount",
        "locusAlleleCount",
        "alleleFrequency",
        "ancillaryResults",
    ]
    for field in expected_fields:
        assert field in caf, f"expected field {field} in CAF"


def test_caf_gen_one_pheno(vrs_id, vcf_path, phenotype_table):
    """test caf generation with default parameters (all patients, all phenotypes)"""

    phenotype = "HP:0001822"
    caf = get_cohort_allele_frequency(
        vrs_id, vcf_path, phenotype_table=phenotype_table, phenotype=phenotype
    )
    print(json.dumps(caf))

    expected_allele_freq = 0.0015
    actual_allele_freq = round(caf["alleleFrequency"], 4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"


# def test_caf_gen_multi_participant():

# def test_caf_gen_one_pheno_multi_participant():

# def test_caf_gen_one_pheno_multi_participant():
