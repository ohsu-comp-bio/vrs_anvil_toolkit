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


@pytest.fixture
def remote_vcf_path():
    return "gs://fc-secure-0392c60d-b346-4e60-b5f1-602a9bbfb0e1/gregor_joint_callset_022824/chr_vcfs/gregor_consortium_dataset_u03_v2.chrY.vcf.gz"


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
        "GCS_OAUTH_TOKEN is required "
        "export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token) "
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


def test_caf_generation_from_vcf_row(vcf_path):
    vcf = VariantFile(vcf_path)
    ref_computed = True
    caf_dicts = []

    for i, record in enumerate(vcf.fetch()):
        vrs_ids = record.info["VRS_Allele_IDs"]
        print("~~~~~~~~~~", i, "~~~~~~~~~~")
        print("....VRS IDs:", vrs_ids)

        focus_allele_count = 0
        locus_allele_count = 0
        num_homozygotes = 0
        num_hemizygotes = 0

        assert len(vcf.header.samples) == len(
            record.samples
        ), "row has different number of samples as the record"

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


###########
# METHODS #
###########


# FIXME: move out of tests
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


def get_patient_phenotype_index(
    phenotypes_table: str = None, cached_dict: str = None, as_set: bool = False
) -> dict[str, list | set]:
    """get list of phenotypes associated with a patient"""

    # load from cache if exists
    if cached_dict is not None and os.path.exists(cached_dict):
        with open(cached_dict, "r") as file:
            patient_to_phenotypes = json.load(file)
            return patient_to_phenotypes

    # if unspecified, parse table using Terra table
    if phenotypes_table is None:
        assert False, "pulling from Terra phenotypes table not implemented yet"

    # table path specified, parse using that table
    phenotype_index = {}
    with open(phenotypes_table, "r") as file:
        if phenotypes_table.endswith(".csv"):
            pheno_df = pd.read_csv(file)
        elif phenotypes_table.endswith(".tsv"):
            pheno_df = pd.read_csv(file, sep="\t")
        else:
            raise Exception(
                "Only csv and tsv file types implemented for phenotype table"
            )

        for participant_id in pheno_df["participant_id"].unique():
            phenotypes = pheno_df[pheno_df["participant_id"] == participant_id][
                "term_id"
            ].unique()

            phenotype_index[participant_id] = (
                set(phenotypes) if as_set else list(phenotypes)
            )

        return phenotype_index


def get_vcf_row(variant_id: str, vcf: VariantFile) -> VariantRecord:
    # TODO: implement, could have multiple VRS IDs to look fofr

    assert "VRS_Allele_IDs" in vcf.header.info, (
        "no VRS_Allele_IDs key in INFO found, "
        "please ensure that this is an VRS annotated VCF"
    )

    for i, record in enumerate(vcf.fetch()):
        if variant_id in record.info["VRS_Allele_IDs"]:
            print(";;;;;;;; record found at", i)
            return record


def get_cohort_allele_frequency(
    variant_id: str,
    vcf_path: str,
    phenotype_table: str = None,
    participant_list: list[str] = None,
    phenotype: str = None,
) -> dict:
    """Get CAF evidence objects given a variant (VRS ID) and optionally
    a list of patients to subset with (defaults to all patients)
    and a list of phenotypes to subset with (defaults to all phenotypes)"""

    assert (
        "ga4gh:VA" in variant_id
    ), "variant ID type not yet supported, use VRS ID instead"

    # get index of variant to patient
    # in this case, the VCF row of the variant_id
    vcf = VariantFile(vcf_path)
    record = get_vcf_row(variant_id, vcf)

    # get index of participant to phenotypes
    print(phenotype_table)
    phenotype_index = get_patient_phenotype_index(phenotype_table, as_set=True)

    # create subsets
    participant_set = (
        set(participant_list) if participant_list is not None else set(record.samples)
    )

    # variables for final cohort allele frequency (CAF) object
    focus_allele_count = 0
    locus_allele_count = 0
    num_homozygotes = 0
    num_hemizygotes = 0
    cohort_phenotypes = set()

    # aggregate data for CAF if...
    for participant_id, genotype in record.samples.items():
        # 1. patient has allele
        alleles = genotype.allele_indices
        if alleles == (None, None):
            continue

        # 2. patient in subset
        if participant_id not in participant_set:
            continue

        # 3. phenotype in phenotypes subset
        if participant_id not in phenotype_index:
            # no matching phenotype for the participant should factor in
            # as 0 of 2 alleles present
            locus_allele_count += 2
            continue

        phenotype_set = phenotype_index[participant_id]
        if phenotype is not None and phenotype not in phenotype_set:
            locus_allele_count += 2
            continue

        # increment allele counts
        num_alleles = sum(alleles)
        focus_allele_count += num_alleles
        locus_allele_count += 2

        # record zygosity
        if num_alleles == 1:
            num_hemizygotes += 1
        elif num_alleles == 2:
            num_homozygotes += 1

        # add to set of all phenotypes
        cohort_phenotypes.update(phenotype_set)

    # populate finaly caf dict
    allele_frequency = focus_allele_count * 1.0 / locus_allele_count

    caf_dict = {
        "type": "CohortAlleleFrequency",
        "label": f"Overall Cohort Allele Frequency for {variant_id}",
        "derivedFrom": {
            "id": "GREGOR_COMBINED_CONSORTIUM_U07",
            "type": "DataSet",
            "label": "GREGoR Combined Consortium U07",
            "version": f"Created {datetime.now()}",
        },
        "focusAllele": variant_id,
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
            "phenotypes": list(cohort_phenotypes),
        },
    }

    return caf_dict


@pytest.fixture()
def vrs_id(vcf_path):
    for i, record in enumerate(VariantFile(vcf_path)):
        if i == 3:
            return record.info["VRS_Allele_IDs"][-2]


@pytest.fixture()
def phenotype_table():
    assert (
        "PHENOTYPE_TABLE" in os.environ
    ), "no PHENOTYPE_TABLE bash variable defining the path of the phenotype table, make sure to export"

    return os.environ["PHENOTYPE_TABLE"]


def test_default_caf_gen(vrs_id, vcf_path, phenotype_table):
    """test caf generation with default parameters (all patients, all phenotypes)"""
    caf = get_cohort_allele_frequency(vrs_id, vcf_path, phenotype_table=phenotype_table)
    print(json.dumps(caf))

    assert round(caf["alleleFrequency"], 6) == 0.040119, "incorrect allele frequency"
