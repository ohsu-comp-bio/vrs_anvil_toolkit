import io
import json
import os
import pandas as pd

from datetime import datetime
from firecloud import api as fapi
from pysam import VariantFile, VariantRecord


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


def get_patient_phenotype_index(
    phenotype_table: str = None, cached_dict: str = None, as_set: bool = False
) -> dict[str, list | set]:
    """get list of phenotypes associated with a patient""

    Args:
        patient_id (str): patient identifier used by phenotypes table
        phenotype_table (str, optional): Path to csv/tsv of phenotype data specified by the GREGoR data model.
                Defaults to pulling from a Terra data table in existing workspace titled "phenotypes".
                For more info on the data model, see https://gregorconsortium.org/data-model
        cached_dict (str, optional): Path to cached dictionary to use. Defaults to None.

    Returns:
        list: list of phenotype codes
    """

    # load from cache if exists
    if cached_dict is not None and os.path.exists(cached_dict):
        with open(cached_dict, "r") as file:
            patient_to_phenotypes = json.load(file)
            return patient_to_phenotypes

    if phenotype_table is None:

        # if unspecified, ensure valid Terra environment
        for env_key in ["WORKSPACE_NAMESPACE", "WORKSPACE_NAME"]:
            assert (
                env_key in os.environ
            ), f"ERROR: No {env_key} key found in environmental variables. Make sure these are set in Terra"

        # create dataframe from Terra data table
        # https://github.com/broadinstitute/fiss/blob/master/firecloud/api.py
        try:
            response = fapi.get_entities_tsv(
                os.environ["WORKSPACE_NAMESPACE"],
                os.environ["WORKSPACE_NAME"],
                "phenotype",
                model="flexible",
            )
            response.raise_for_status()
        except Exception as e:
            if response.json() and e in response.json():
                error_message = response.json()["message"]
            else:
                error_message = e
            print(
                f"Error while loading phenotype data table from workspace: \n{error_message}"
            )

        phenotype_tsv = io.StringIO(response.text)
        phenotype_df = pd.read_csv(phenotype_tsv, sep="\t")
    else:
        # table path specified, parse using that table
        phenotype_index = {}
        with open(phenotype_table, "r") as file:
            if phenotype_table.endswith(".csv"):
                phenotype_df = pd.read_csv(file)
            elif phenotype_table.endswith(".tsv"):
                phenotype_df = pd.read_csv(file, sep="\t")
            else:
                raise Exception(
                    "Only csv and tsv file types implemented for phenotype table"
                )

    for participant_id in phenotype_df["participant_id"].unique():
        phenotypes = phenotype_df[phenotype_df["participant_id"] == participant_id][
            "term_id"
        ].unique()

        phenotype_index[participant_id] = (
            set(phenotypes) if as_set else list(phenotypes)
        )

    return phenotype_index


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

    # get index of participant to phenotypes
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
