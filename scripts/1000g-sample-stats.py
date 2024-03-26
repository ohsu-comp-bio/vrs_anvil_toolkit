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
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import os
import pandas as pd
import pysam
import seaborn as sns
import yaml

from collections import defaultdict
from vrs_anvil import query_metakb


# settings and files to load
# yaml_path = "state/metrics_20240320_143108.yaml" # chr1 results from 1000g chr1...vcf.gz
# yaml_path = "state/metrics_20240321_095608.yaml"  # chr1 and 2
yaml_path = "/Users/wongq/Downloads/metrics_20240326_160440.yaml"
figure_dir = "figures"
variant_percentages_file_name = "chr1_to_chr2.png"
variant_histogram_file_name = "variants_per_patient"
save_figures = False
show_figures = True

# seaborn styling
sns.set_theme()
sns.set_style("whitegrid")

with open(yaml_path, "r") as file:
    metrics = yaml.safe_load(file)

# data storage
sample_dict = defaultdict(list)

# get evidence key drill down for each ting
METAKB_DIR = f"../tests/fixtures/metakb"
json_paths = list(pathlib.Path(METAKB_DIR).glob("*.json"))

for file_path in metrics:
    if file_path == "total":
        continue

    # load original vcf
    real_path = os.path.realpath(file_path)
    vcf_reader = pysam.VariantFile(real_path)

    # check if metakb evidence found
    print("from file", file_path)
    assert "evidence" in metrics[file_path], "no evidence found"
    evidence = metrics[file_path]["evidence"]

    # collect info for each vrs allele id
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
                if any(
                    int(allele) > 0 for allele in gt_values
                ):  # Check if any allele is non-zero
                    id_count += 1
                    if sample not in sample_dict:
                        sample_dict[sample] = {}

                    if "vrs_ids" in sample_dict[sample]:
                        sample_dict[sample]["vrs_ids"].append(allele_id)
                        sample_dict[sample]["study_ids"].extend(study_ids)
                    else:
                        sample_dict[sample]["vrs_ids"] = [allele_id]
                        sample_dict[sample]["study_ids"] = list(study_ids)

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
metakb_percent = get_percent(len(sample_dict), num_samples, output=True)

knowledgebases = ["MOAlmanac", "CIVIC", "All Knowledgebases"]
keywords = ["moa", "civic"]
all_percentages = []

for keyword in keywords:
    count = sum(
        [
            any(keyword in study_id for study_id in id_lists["study_ids"])
            for _, id_lists in sample_dict.items()
        ]
    )
    all_percentages.append(count * 100.0 / num_samples)


all_percentages.append(metakb_percent)

# bar plot
sns.barplot(x=knowledgebases, y=all_percentages, hue=knowledgebases)
sns.despine(left=True)

# Set labels and title
plt.xlabel("Knowledgebases")
plt.ylabel("Percent of Patients with Match (%)")
plt.title(
    "[chr1 to chr2] Percent of 1000G Patients with a Variant Match to a Knowledgebase",
    wrap=True,
)
plt.ylim(0, 100)
plt.yticks(range(0, 101, 20))

# Show the plot
plt.tight_layout()
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
