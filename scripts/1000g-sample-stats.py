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
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import os
import pysam
import vcf
import yaml

from collections import defaultdict
from glom import glom


# yaml_path = "state/metrics_20240320_143108.yaml" # chr1 results from 1000g chr1...vcf.gz
yaml_path = "state/metrics_20240321_095608.yaml"  # chr1 and 2
figure_dir = "figures"
figure_name = "chr1_to_chr2.png"
save_figure = True
show_figure = True

with open(yaml_path, "r") as file:
    metrics = yaml.safe_load(file)


def find_keys_from_ids(data, target="ga4gh:VA", prefix="", result=None):
    if result is None:
        result = set()

    if isinstance(data, dict):
        for key, value in data.items():
            new_prefix = f"{prefix}.{key}" if prefix else key
            if isinstance(value, (dict, list)):
                find_keys_from_ids(value, target, new_prefix, result)
            elif isinstance(value, str) and target in value:
                result.add(new_prefix)
    elif isinstance(data, list):
        for index, value in enumerate(data):
            new_prefix = f"{prefix}.{index}"
            if isinstance(value, (dict, list)):
                find_keys_from_ids(value, target, new_prefix, result)
            elif isinstance(value, str) and target in value:
                result.add(new_prefix)
    return result


# data storage
sample_dict = defaultdict(list)
allele_id_to_count = defaultdict(int)

# get evidence key drill down for each ting
METAKB_DIR = f"../tests/fixtures/metakb"
json_paths = list(pathlib.Path(METAKB_DIR).glob("*.json"))

for file_path in metrics:
    if file_path == "total":
        continue

    # load original vcf
    real_path = os.path.realpath(file_path)
    vcf_reader = pysam.VariantFile(real_path)

    print("from file", file_path)
    if "evidence" not in metrics[file_path]:
        print("no evidence found")
        continue
    evidence = metrics[file_path]["evidence"]

    for allele_id, allele_info in evidence.items():
        gnomad_expr = allele_info["parameters"]["var"]
        chrom, pos, ref, alt = gnomad_expr.split("-")
        pos = int(pos)

        # should be a single record each time
        assert (
            sum([1 for _ in vcf_reader.fetch("chr1", pos - 1, pos)]) == 1
        ), f"more than one record found at pos {pos}"

        for record_idx, record in enumerate(vcf_reader.fetch("chr1", pos - 1, pos)):
            print(f"gnomad: {gnomad_expr}")
            print(f"allele_id: {allele_id}")

            assert ref == record.ref, f"expected ref {record.ref}, got {ref}"
            assert alt in record.alts, f"expected one of {record.alt} alts, got {alt}"

            # get list of knowledgebase keys (eg civic.mpid:254)
            # TODO: should this be a part of the cache?
            study_ids = []

            def get_values_from_keys(data, keys, final_key):
                final_suffix = "." + final_key
                return [
                    glom(data, ".".join(key.split(".")) + final_suffix) for key in keys
                ]

            # get study id associated with vrs allele id
            for json_path in json_paths:
                with open(json_path) as file:
                    data = json.load(file)

                allele_id_keys = find_keys_from_ids(data, allele_id)

                json_study_keys = [
                    ".".join(key.split(".")[:2]) for key in allele_id_keys
                ]

                json_study_ids = get_values_from_keys(data, json_study_keys, "id")
                study_ids.extend(json_study_ids)

                # print id and descriptions
                for study_key, study_id in zip(json_study_keys, json_study_ids):
                    if study_key.startswith("studies"):
                        final_suffix = ".description"
                    elif study_key.startswith("molec"):
                        final_suffix = ".aliases"

                    description = glom(
                        data, ".".join(study_key.split(".")) + final_suffix
                    )
                    print(f"{study_id}: {description}")

            # get related stuff
            for sample, genotype in record.samples.items():
                gt_values = genotype["GT"]
                if any(
                    int(allele) > 0 for allele in gt_values
                ):  # Check if any allele is non-zero
                    allele_id_to_count[allele_id] += 1
                    if sample not in sample_dict:
                        sample_dict[sample] = {}

                    if "vrs_ids" in sample_dict[sample]:
                        sample_dict[sample]["vrs_ids"].append(allele_id)
                        sample_dict[sample]["study_ids"].extend(study_ids)
                    else:
                        sample_dict[sample]["vrs_ids"] = [allele_id]
                        sample_dict[sample]["study_ids"] = list(study_ids)

            print(f"total count: {allele_id_to_count[allele_id]}\n")

    # use glom to get the key

#### Figures ####
# average number of variants per individual


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
colors = ["skyblue", "lightcoral", "gold", "grey"]

# bar plot
plt.bar(knowledgebases, all_percentages, color=mcolors.TABLEAU_COLORS)

# Set labels and title
plt.xlabel("Knowledgebases")
plt.ylabel("Percent of Patients with Match")
plt.title(
    "[chr1 to chr2] Percent of 1000G Patients with a Variant Match to a Knowledgebase",
    wrap=True,
)
plt.ylim(0, 120)
plt.yticks(range(0, 110, 20))

# Show the plot
# plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
os.makedirs(figure_dir, exist_ok=True)
if save_figure:
    plt.savefig(f"{figure_dir}/{figure_name}", dpi=300)
if show_figure:
    plt.show()


# # load the line number and
# def get_record_at_line(vcf_filename, pos):
#     # Open the VCF file with indexing enabled
#     vcf_file = pysam.VariantFile(vcf_filename, 'r', index_filename=vcf_filename + '.tbi')

#     # Retrieve the record at the specified line number
#     try:
#         record = vcf_file.fetch(start=pos, end=pos)
#         return next(record)  # Return the first record found
#     except StopIteration:
#         return None  # Line number not found
