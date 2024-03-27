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
from vrs_anvil import query_metakb, truncate


# settings and files to load
# yaml_path = "state/metrics_20240320_143108.yaml" # chr1 results from 1000g chr1...vcf.gz
metrics_path = "state/metrics_20240321_095608.yaml"  # chr1 and 2
figure_dir = "figures"
variant_percentages_file_name = "chr1_to_chr2.png"
variant_histogram_file_name = "chr1_to_chr2_variants_per_patient"
save_figures = True
show_figures = False

# seaborn styling
sns.set_theme()
sns.set_style("whitegrid")

# data storage
sample_dict = defaultdict(list)
variant_types_set = set()

# get evidence key drill down for each ting
METAKB_DIR = f"../tests/fixtures/metakb"
json_paths = list(pathlib.Path(METAKB_DIR).glob("*.json"))

evidence_per_file = {}

# collect all evidence from a list of metrics files
# chrs = [str(num) for num in range(1,23)] + ["X"]
chrs = [str(num) for num in range(1,3)]
for c in chrs:
    all_metrics_paths = glob(f"all-chrs/chr{c}/state/metrics*.yaml")
    most_recent_metrics_path = all_metrics_paths[-1]

    with open(most_recent_metrics_path, "r") as file:
        metrics = yaml.safe_load(file)

    for file_path in metrics:
        if file_path == "total":
            continue

        # load original vcf
        original_path = os.path.realpath(file_path)
        vcf_reader = pysam.VariantFile(original_path)

        # check if metakb evidence found
        print("from file", truncate(original_path, 0, 46))
        if "evidence" in metrics[file_path]:
            print(f"{len(metrics[file_path]["evidence"])} hits found")
            evidence_per_file.update({original_path: metrics[file_path]["evidence"]})
        else:
            print("no evidence from this file")
        print()

# collect info for each vrs allele id
for file_path, evidence in evidence_per_file.items():
    vcf_reader = pysam.VariantFile(file_path)

    for allele_id, allele_info in evidence.items():

        # extract coordinate information
        gnomad_expr = allele_info["parameters"]["var"]
        chrom, pos, ref, alt = gnomad_expr.split("-")
        pos = int(pos)
        print(f"\tgnomad: {gnomad_expr}")
        print(f"\tallele_id: {allele_id}")

        # get study id associated with vrs allele id
        metakb_response = query_metakb(allele_id, log=True)
        assert metakb_response is not None, f"no metakb hit for allele {allele_id}"
        study_ids = metakb_response["study_ids"]

        assert study_ids == [study["id"] for study in metakb_response["studies"]], \
            "study ids in wrong order to add variant types"
        variant_types = [study["qualifiers"]["alleleOrigin"] for study in metakb_response["studies"]]
        variant_types_set.update(variant_types)

        for study in metakb_response["studies"]:
            print(f"\t{study['type']} ({study['id']}): {study['description']}")

        # should be a single record each time
        assert (
            sum([1 for _ in vcf_reader.fetch(chrom, pos - 1, pos)]) == 1
        ), f"more than one record found at pos {pos}"

        for record_idx, record in enumerate(vcf_reader.fetch(chrom, pos - 1, pos)):

            assert ref == record.ref, f"expected ref {record.ref}, got {ref}"
            assert alt in record.alts, f"expected one of {record.alts} alts, got {alt}"

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
        print(type(percent))
        if data is None:
            data = np.array([knowledgebases[i], percent, variant_type])
        else:
            data = np.vstack([data, [knowledgebases[i], percent, variant_type]])

expected_num_rows = len(knowledgebases) * len(variants)
assert len(data) == expected_num_rows, f"expected {expected_num_rows} rows, got {len(data)}"

# bar plot
df_pct = pd.DataFrame(data, columns=cols).astype(dtype={col: dtype for col, dtype in zip(cols, dtypes)})
print(df_pct)
print(df_pct[cols[PCT]])
print(df_pct.dtypes)
ax = sns.barplot(data=df_pct, x=cols[VAR], y=cols[PCT], hue=cols[KB], order=["somatic", "germline", TOTAL])

# # add percent labels
# for i, bar in enumerate(ax.patches):
#     height = bar.get_height()
#     label_text = f"{df_pct.iloc[cols[KB], i]:.1f}%" # Format as percentage with one decimal place
#     label_x = bar.get_x() + bar.get_width() / 2  # Center the label horizontally
#     label_y = bar.get_y() + height / 2         # Position slightly above the bar
#     ax.text(label_x, label_y, label_text, ha='center', va='center')

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
plt.legend(loc='upper left')


# Show the plot
# plt.tight_layout()
os.makedirs(figure_dir, exist_ok=True)
if save_figures:
    plt.savefig(f"{figure_dir}/{variant_percentages_file_name}", dpi=300)
if show_figures:
    plt.show()

# average number of variants for all samples
num_variants_per_patient = [len(v["vrs_ids"]) for v in sample_dict.values()]
num_variants_per_patient.extend([0 for _ in range(num_samples - len(sample_dict))])

NUM_VARIANTS = "num_variants"

plt.figure()
df = pd.DataFrame({NUM_VARIANTS: num_variants_per_patient})
df["percentage"] = df[NUM_VARIANTS].value_counts(normalize=True) * 100
sns.histplot(data=df, x=NUM_VARIANTS, stat="density", bins=3, discrete=True)
plt.xlabel("Number of Variants")
plt.ylabel("Percentage of All Patients (%)")
plt.title("Number of Variants Associated with Each Patient")

plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(1))
plt.grid(False)
plt.xticks(range(min(df[NUM_VARIANTS]), max((df[NUM_VARIANTS]))))

if save_figures:
    plt.savefig(f"{figure_dir}/{variant_histogram_file_name}", dpi=300)
if show_figures:
    plt.show()
