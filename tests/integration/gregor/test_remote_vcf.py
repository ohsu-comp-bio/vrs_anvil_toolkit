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


def test_caf_generation_from_vcf_row(vcf_path):
    vcf = VariantFile(vcf_path)
    ref_computed = True
    caf_dicts = []

    assert "VRS_Allele_IDs" in vcf.header.info, (
        "no VRS_Allele_IDs key in INFO found, "
        "please ensure that this is an VRS annotated VCF"
    )

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
    assert "VRS_Allele_IDs" in vcf.header.info, (
        "no VRS_Allele_IDs key in INFO found, "
        "please ensure that this is an VRS annotated VCF"
    )

    # TODO: use VRS ID -> VCF index to grab
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
    """Create a cohort allele frequency for either genotypes or phenotypes

    Args:
        variant_id (str): variant ID (VRS ID)
        vcf_path (str): path to VCF
        phenotype_table (str, optional): where to pull phenotype information from. Defaults to None.
        participant_list (list[str], optional): Subset of participants to use. Defaults to None.
        phenotype (str, optional): Specific phenotype to subset on. Defaults to None.

    Returns:
        dict: Cohort Allele Frequency object
    """

    assert (
        "ga4gh:VA" in variant_id
    ), "variant ID type not yet supported, use VRS ID instead"

    # get index of variant to patient
    # in this case, the VCF row of the variant_id
    vcf = VariantFile(vcf_path)
    record = get_vcf_row(variant_id, vcf)

    # if multiple alts, get index associated with alt
    # if ref is specified in VRS Alleles IDs, adjust indices to match
    alt_index = record.info["VRS_Allele_IDs"].index(variant_id)
    if "REF" not in vcf.header.info["VRS_Allele_IDs"].description:
        alt_index -= 1
    print("alt_index:", alt_index)

    # get index of participant to phenotypes
    print(phenotype_table)
    phenotype_index = get_patient_phenotype_index(phenotype_table, as_set=True)

    # create cohort, defaults to all patients in VCF
    cohort = (
        set(participant_list) if participant_list is not None else set(record.samples)
    )

    # FIXME: support for hemizygous regions (chrX / mitochondrial variants)
    # GREGOR makes use of DRAGEN's continuous allele frequency approach:
    # https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/MitochondrialCalling.htm
    within_hemizygous_region = record.chrom in ["chrM", "chrY"]
    within_x_chr = record.chrom == "chrX"

    # variables for final cohort allele frequency (CAF) object
    focus_allele_count = 0
    locus_allele_count = 0
    cohort_phenotypes = set() if phenotype is None else [phenotype]

    # if haploid not diploid, these fields are not relevant
    if len(record.samples[0]["GT"]) == 2:
        num_homozygotes = 0
        num_hemizygotes = 0

    # aggregate data for CAF so long as...
    for participant_id, genotype in record.samples.items():
        # 1. participant has recorded genotype
        alleles = genotype.allele_indices
        if all(a is None for a in alleles):
            continue

        # 2. participant in specified cohort
        if participant_id not in cohort:
            continue

        # 3. participant did not specify phenotype (doing general query)
        # or participant has specified phenotype
        has_specified_phenotype = (
            participant_id in phenotype_index
            and phenotype in phenotype_index[participant_id]
        )

        # with these conditions satisfied...
        if phenotype is None or has_specified_phenotype:
            # increment focus allele count, handling multiple alts edge case
            num_alt_alleles = sum(
                [1 for _, alt_number in enumerate(alleles) if alt_number == alt_index]
            )

            if not within_hemizygous_region:
                focus_allele_count += num_alt_alleles
            elif num_alt_alleles > 0:
                focus_allele_count += 1

            # record zygosity
            if not within_hemizygous_region and len(alleles) == 2:
                if num_alt_alleles == 1:
                    num_hemizygotes += 1
                elif num_alt_alleles == 2:
                    num_homozygotes += 1

        # increment total allele count
        if within_hemizygous_region:
            locus_allele_count += 1
        elif within_x_chr:
            # FIXME: make use of sex of participant?
            is_female = True
            locus_allele_count += 2 if is_female else 1
        else:
            locus_allele_count += len(alleles)

        # update phenotypes as necessary
        if phenotype is not None:
            # # add to set of all phenotypes
            # if participant_id in phenotype_index:
            #     cohort_phenotypes.update(phenotype_index[participant_id])
            continue
        else:
            # aggregate phenotypes if they exist
            if participant_id in phenotype_index:
                cohort_phenotypes.update(phenotype_index[participant_id])

    # populate final caf dict
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
            return record.info["VRS_Allele_IDs"][1]  # 1 since 0 is a ref


@pytest.fixture()
def phenotype_table():
    assert (
        "PHENOTYPE_TABLE" in os.environ
    ), "no PHENOTYPE_TABLE bash variable defining the path of the phenotype table, make sure to export"

    return os.environ["PHENOTYPE_TABLE"]


def test_allele_counts_first_10_rows(phenotype_table):
    vcf_path = "/Users/wongq/data/gregor/ann.100k.chr3.vcf.gz"
    vcf = VariantFile(vcf_path)
    has_ref = "REF" in vcf.header.info["VRS_Allele_IDs"].description

    for i, record in enumerate(vcf.fetch()):
        print("~~~~~~~~~ row ", i, "~~~~~~~~~")

        # if REF VRS ID recorded, ignore
        vrs_allele_ids = record.info["VRS_Allele_IDs"]
        print("vrs_allele_ids:", vrs_allele_ids)
        if has_ref:
            vrs_allele_ids = vrs_allele_ids[1:]

        for alt_index, allele_id in enumerate(vrs_allele_ids):
            print("alt id:", allele_id)
            if not allele_id:
                continue

            caf = get_cohort_allele_frequency(
                allele_id, vcf_path, phenotype_table=phenotype_table
            )

            print("alt_index", alt_index)
            print("AC:", record.info["AC"], caf["focusAlleleCount"])
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
    phenotype = "HP:0001822"
    """test caf generation with default parameters (all patients, all phenotypes)"""
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
