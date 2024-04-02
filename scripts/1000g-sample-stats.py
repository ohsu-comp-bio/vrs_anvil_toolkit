"""
Get statistics from 1000g samples, including
  - total samples
  - % matches to each kb individually and metakb as a whole
  - study ids and descriptions by variant hit in metakb

caveats: variant match on vrs id, not on
  - matching zygosity
  - germline vs somatic variant

(to be run in a tmp/ directory)
"""

import json
import pathlib
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
import yaml

from collections import defaultdict
from glob import glob
from vrs_anvil import query_metakb

# grabs all metrics files from this directory
metrics_dir = "state/chr1_and_2"
# metrics_path = "state/chr1" # chr1 results from 1000g chr1...vcf.gz

# settings and files to load
figure_dir = "figures"
variant_percentages_file_name = "chr1_to_chr2.png"
variant_histogram_file_name = "chr1_to_chr2_variants_per_patient"
save_figures = False
show_figures = True

# cohort allele frequency info
caf_label_name = "vrs-python"
caf_label_version = "2.0.0-a5"
caf_dicts = []

# seaborn styling
sns.set_theme()
sns.set_style("whitegrid")

# helper function
def truncate(s, first_few, last_few):
    "truncate string printing only first_few and last_few characters"
    return f"{s[:first_few]}...{s[-last_few:]}"


# data storage
sample_dict = defaultdict(list)
variant_types_set = set()

# get evidence key drill down for each ting
METAKB_DIR = f"../tests/fixtures/metakb"
json_paths = list(pathlib.Path(METAKB_DIR).glob("*.json"))

evidence_per_file = {}

# collect all evidence from a list of metrics files
# chrs = [str(num) for num in range(1,23)] + ["X"]
chrs = [str(num) for num in range(1, 3)]
for c in chrs:
    all_metrics_paths = glob(f"{metrics_dir}/metrics*.yaml")
    most_recent_metrics_path = all_metrics_paths[-1]

    with open(most_recent_metrics_path, "r") as file:
        metrics = yaml.safe_load(file)

    for file_path in metrics:
        if file_path == "total":
            continue

        # load original vcf
        original_path = os.path.realpath(file_path)
        vcf_reader = pysam.VariantFile(original_path)
        num_samples = len(vcf_reader.header.samples)
        assert num_samples == 3202, \
            f"expected 3202 samples for {file_path}, got {num_samples}"

        # check if metakb evidence found
        print("from file", truncate(original_path, 0, 46))
        if "evidence" in metrics[file_path]:
            print(f"{len(metrics[file_path]['evidence'])} hits found")
            evidence_per_file.update({original_path: metrics[file_path]["evidence"]})
        else:
            print("no evidence from this file")
        print()

# get study info (id: description) for each variant
for file_path, evidence in evidence_per_file.items():
    vcf_reader = pysam.VariantFile(file_path)

    print(truncate(file_path, 0, 47))

    for allele_id, allele_info in evidence.items():

        # extract coordinate information
        gnomad_expr = allele_info["parameters"]["var"]
        chrom, pos, ref, alt = gnomad_expr.split("-")
        pos = int(pos)
        print(f"\tgnomad: {gnomad_expr}")
        print(f"\tallele_id: {allele_id}")

        # get study id associated with vrs allele id
        metakb_response = query_metakb(allele_id, log=True)
        if metakb_response is None:
            print(f"no metakb hit for allele {allele_id}\n")
            continue

        study_ids = metakb_response["study_ids"]

        assert study_ids == [
            study["id"] for study in metakb_response["studies"]
        ], "study ids in wrong order to add variant types"
        variant_types = [
            study["qualifiers"]["alleleOrigin"] for study in metakb_response["studies"]
        ]
        variant_types_set.update(variant_types)

        for study in metakb_response["studies"]:
            print(f"\t{study['type']} ({study['id']}): {study['description']}")

        # should be a single record each time
        if sum([1 for _ in vcf_reader.fetch(chrom, pos - 1, pos)]) != 1:
            print(f"[WARNING] more than one record found at pos {pos}")

        for record_idx, record in enumerate(vcf_reader.fetch(chrom, pos - 1, pos)):

            assert ref == record.ref, f"expected ref {record.ref}, got {ref}"
            if alt not in record.alts:
                print(
                    f"Skipping {gnomad_expr}, expected one of {record.alts} alts, got {alt}"
                )

            # populate patient/sample dictionary
            id_count = 0
            for sample, genotype in record.samples.items():
                gt_values = genotype["GT"]

                # Check if any allele is non-zero
                if any(int(allele) > 0 for allele in gt_values):
                    id_count += 1
                    if sample not in sample_dict:
                        sample_dict[sample] = {}

                    if "vrs_ids" in sample_dict[sample]:
                        sample_dict[sample]["vrs_ids"].append(allele_id)
                        sample_dict[sample]["study_ids"].extend(study_ids)
                        sample_dict[sample]["variant_types"].extend(variant_types)
                    else:
                        sample_dict[sample]["vrs_ids"] = [allele_id]
                        sample_dict[sample]["study_ids"] = list(study_ids)
                        sample_dict[sample]["variant_types"] = list(variant_types)

            print(f"\ttotal count: {id_count}\n")

            # Create cohort allele frequency objects
        # get counts
        num_patients_w_allele = id_count
        allele_freq = num_patients_w_allele*1.0/num_samples

        # create caf dict
        # TODO: locusAlleleCount: The total number of alleles observed at
        # that specific gene locus in the cohort (including the focus allele
        # and any alternate alleles)

        caf_dict = {
            "id": allele_id,
            "type": "CohortAlleleFrequency",
            "label": f"Cohort Allele Frequency for {allele_id} ({gnomad_expr})",
            "derivedFrom": {
                "id": f"AnVIL_1000G_PRIMED-data-model",
                "type": "DataSet",
                "label": f"AnVIL 1000G PRIMED Data Model",
                "version": "Last Updated 3/12/2024",
            },
            "focusAllele": {
                "oneOf": [allele_id],
            },
            "focusAlleleCount": num_patients_w_allele,
            "locusAlleleCount": num_samples,
            "alleleFrequency": allele_freq,
            "cohort": {
                "id": "all_3202_samples",
                "label": "High-coverage sequencing of 3202 samples",
            }
        }
        caf_dicts.append(caf_dict)


print("Cohort Allele Frequency Objects:")
for caf_dict in caf_dicts:
    print(json.dumps(caf_dict, indent=2), "\n")


#### Figures ####
# patients with one match total
def get_percent(a, b, output=True):
    "pretty print percentages"
    percent = float(f"{(100.0*a/b):.1f}")
    if output:
        print(f"{a}/{b} = {percent}%")
    return percent


print("patients with at least one variant match:", end=" ")
num_samples = len(vcf_reader.header.samples)

knowledgebases = ["MOAlmanac", "CIVIC", "All Knowledgebases"]
kb_keywords = ["moa", "civic", ""]  # "" represents match for everything

KB, PCT, VAR = range(3)
cols = ["knowledgebase", "percent", "variant_type"]
dtypes = [str, float, str]

TOTAL = "all"
variants = sorted(list(variant_types_set)) + [TOTAL]
data = None

# for each variant type (germline, somatic) and kb, get count of people in
for variant_type in variants:
    for i, keyword in enumerate(kb_keywords):
        num_matching_samples = 0

        # increment if sample has matching variant and kb
        for _, id_lists in sample_dict.items():
            v_types = np.array(id_lists["variant_types"])

            if variant_type == TOTAL:
                study_ids = np.array(id_lists["study_ids"])
            else:
                study_ids = np.array(id_lists["study_ids"])[v_types == variant_type]
            if any(keyword in study_id for study_id in study_ids):
                num_matching_samples += 1

        percent = num_matching_samples * 100.0 / num_samples
        if data is None:
            data = np.array([knowledgebases[i], percent, variant_type])
        else:
            data = np.vstack([data, [knowledgebases[i], percent, variant_type]])

expected_num_rows = len(knowledgebases) * len(variants)
assert (
    len(data) == expected_num_rows
), f"expected {expected_num_rows} rows, got {len(data)}"

# bar plot
df_pct = pd.DataFrame(data, columns=cols).astype(
    dtype={col: dtype for col, dtype in zip(cols, dtypes)}
)
ax = sns.barplot(
    data=df_pct,
    x=cols[VAR],
    y=cols[PCT],
    hue=cols[KB],
    order=["somatic", "germline", TOTAL],
)


# add percent labels
def add_percent_labels(num_bars, ax, vertical_nudge, from_decimal=False):
    multiplier = 100 if from_decimal else 1

    for i, bar in enumerate(ax.patches):
        if i >= num_bars:
            break

        height = bar.get_height()
        percent = f"{multiplier * height:.2f}"
        x = bar.get_x() + bar.get_width() / 2
        y = bar.get_y() + height + vertical_nudge
        ax.text(x, y, percent, ha="center", va="center", fontsize=10)


add_percent_labels(len(df_pct), ax, 2.5)

# sns.despine(left=True)
plt.ylim(0, 100)
plt.yticks(np.arange(0, 101, 20))

# Set labels and title
plt.xlabel("Knowledgebases")
plt.ylabel("Percent of Patients with Match (%)")
plt.title(
    "[chr1 to chr2] Percent of 1000G Patients with a Variant Match (n=3202)",
    wrap=True,
)
plt.legend(loc="upper left")

# Show the plot
plt.tight_layout()
os.makedirs(figure_dir, exist_ok=True)
if save_figures:
    plt.savefig(f"{figure_dir}/{variant_percentages_file_name}", dpi=300)
if show_figures:
    plt.show()

# average number of variants for all samples
num_variant_hits = sum(len(e) for e in evidence_per_file.values())
print("total hits across all patients:", num_variant_hits)

jitter = [0.008, 0.018, 0.004]


NUM_VARIANTS = "num_variants"

for i, variant_type in enumerate(variants):
    num_variants_per_patient = [
        len(
            [
                v
                for v in values["variant_types"]
                if v == variant_type or variant_type == TOTAL
            ]
        )
        for values in sample_dict.values()
    ]

    num_variants_per_patient.extend([0 for _ in range(num_samples - len(sample_dict))])
    assert (
        len(num_variants_per_patient) == num_samples
    ), f"only {len(num_variants_per_patient)} samples, expected {num_samples}"

    df = pd.DataFrame({NUM_VARIANTS: num_variants_per_patient})
    df["percentage"] = df[NUM_VARIANTS].value_counts(normalize=True) * 100

    plt.figure()
    ax = sns.histplot(data=df, x=NUM_VARIANTS, stat="density", discrete=True)
    plt.xlabel("Number of Variants")
    plt.ylabel("Percentage of All Patients (%)")
    plt.title(
        f"Number of {variant_type.capitalize()} Variants Associated with Each Patient"
    )

    variant_hits_by_type = max((df[NUM_VARIANTS]) + 1)
    add_percent_labels(variant_hits_by_type, ax, jitter[i], from_decimal=True)

    plt.gca().yaxis.set_major_formatter(lambda x, _: f"{(x*100):.0f}")
    plt.grid(False)
    plt.xticks(range(min(df[NUM_VARIANTS]), variant_hits_by_type))

    if save_figures:
        plt.savefig(f"{figure_dir}/{variant_histogram_file_name}", dpi=300)
    if show_figures:
        plt.show()
