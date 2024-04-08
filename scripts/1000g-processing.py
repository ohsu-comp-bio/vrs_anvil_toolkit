"""
Get statistics from 1000g samples, including
  - total samples
  - % matches to each kb individually and metakb as a whole
  - study ids and descriptions by variant hit in metakb

(to be run in a tmp/ directory)
"""

import io
import json
import os
import pandas as pd
import pysam
import subprocess
import yaml

from collections import defaultdict
from google.cloud import storage
from glob import glob
from firecloud import api as fapi
from vrs_anvil import query_metakb
from vrs_anvil.annotator import MATCHES, TOTAL, VRS_OBJECT

setup_cmd = f"cd ~/vrs-python-testing && bash terra/setup.sh"
output = subprocess.run(setup_cmd, shell=True, check=True)

INPUT_DIR="/home/jupyter/vrs-python-testing/tmp/work/high_coverage_3202_samples"
SEQREPO_DIR = "/home/jupyter/seqrepo/latest"
BASE_DIR = "/home/jupyter/vrs-python-testing/tmp"

os.makedirs(INPUT_DIR, exist_ok=True)

# helper function
def truncate(s, first_few, last_few):
    "truncate string printing only first_few and last_few characters"
    return f"{s[:first_few]}...{s[-last_few:]}"


# pull down all metrics files stuff

# constants
SOURCE_BUCKET = os.getenv("WORKSPACE_BUCKET").split("//")[1]
TIMESTAMP = "20240329_013953"
BUCKET_DIR = "scattered/state-04-07-2024"

# Find matching files
client = storage.Client()
source_bucket = client.get_bucket(SOURCE_BUCKET)
all_blobs = source_bucket.list_blobs(prefix=BUCKET_DIR)

blobs = [blob for blob in all_blobs if blob.name.endswith(".yaml")]

print("number of metrics files:", len(blobs))


# create output directory
RESULTS_DIR = os.path.expanduser(f"~/vrs-python-testing/tmp/state/{TIMESTAMP}")
os.makedirs(RESULTS_DIR, exist_ok=True)
print(f"downloading to {RESULTS_DIR}:")

# get blob if not already downloaded
metrics_paths = []

for blob in blobs:
    file_name = blob.name.split("/")[-1]
    file_path = f"{RESULTS_DIR}/{file_name}"
    metrics_paths.append(file_path)

    if os.path.exists(file_path):
        print(f"[SKIPPED] {file_name} already exists")
    else:
        blob.download_to_filename(file_path)

        print(f"downloaded {file_path}")

# get original 1000g metadata
df = pd.read_csv(io.StringIO(fapi.get_entities_tsv("anvil-datastorage", \
                "AnVIL_1000G_PRIMED-data-model", "sequencing_file", model="flexible").text), sep='\t')

# get rid of gvcf data
df_vcf = df[df['file_type'].isin(['VCF', 'VCF index'])]
df_1kgp = df_vcf[df_vcf['file_path'].str.contains('1kGP')]

num_vcf_idx_files = sum(df_1kgp['file_type'] == 'VCF index')
num_vcf_files = sum(df_1kgp['file_type'] == 'VCF')
assert num_vcf_files == 23 and num_vcf_idx_files == 23, \
    f"check number of files, {num_vcf_files} vcfs and {num_vcf_idx_files} index files"

# load all vcf files if not pulled
chrs = [str(num) for num in range(1,23)] + ["X"]

if chrs is not None:
    df_chrs = df_1kgp[df_1kgp['chromosome'].isin(set(chrs))]
else:
    df_chrs = df_1kgp
    specific_chromosomes = df_chrs['chromosome'].unique()

assert(len(df_chrs) == 2*len(chrs)), \
    f"Expected 2 files per chr but {len(df_chrs)} files and {len(chrs)} chrs"

# pull from cloud and save paths and chr order
raw_vcfs = []
chrs_of_interest = []

for _, row in df_chrs.iterrows():
    uri = row["file_path"]
    file_name = uri.split("/")[-1]
    output_file = f"{INPUT_DIR}/{file_name}"

    # save vcfs for later, ignore indices
    if row["file_type"] == 'VCF':
        raw_vcfs.append(output_file)
        chrs_of_interest.append(row["chromosome"])

    if os.path.exists(output_file):
        print(f"[SKIPPING] {truncate(file_name, 35, 10)} already exists")
    else:
        pull_vcf_cmd = f"gsutil -u $GOOGLE_PROJECT cp {uri} {INPUT_DIR}/"
        output = subprocess.run(pull_vcf_cmd, shell=True, check=True)

# collect all evidence from a list of metrics files in a directory
matches_per_file = {}
all_metrics_paths = sorted(glob(f"{RESULTS_DIR}/metrics*.yaml"))

total_matches = 0
for metrics_path in all_metrics_paths:

    with open(metrics_path, "r") as file:
        metrics = yaml.safe_load(file)

    # for each file in the metrics, add matches to dict
    for file_path, metrics_dict in metrics.items():
        if file_path == TOTAL:
            continue

        # load original vcf
        original_path = os.path.expanduser(f"{BASE_DIR}/{file_path}")
        vcf_reader = pysam.VariantFile(original_path)
        num_samples = len(vcf_reader.header.samples)
        assert num_samples == 3202, \
            f"expected 3202 samples for {file_path}, got {num_samples}"

        # check if matches to metakb cache found
        print("from file", truncate(original_path, 0, 46))
        if MATCHES in metrics_dict:
            expected_matches = int(metrics_dict["metakb_hits"])
            actual_matches = len(metrics_dict[MATCHES])
            assert expected_matches == actual_matches, \
                f"expected {expected_matches} metakb matches, found {actual_matches}"
            print(f"{actual_matches} matches found")

            total_matches += actual_matches
            matches_per_file.update({original_path: metrics_dict[MATCHES]})
        else:
            print("no matches from this file")
        print()

print(f"total matches: {total_matches}")

# print all files with matches
print("files with evidence:")
for file_path in sorted(matches_per_file.keys()):
    print(truncate(file_path, 0, 47))

# function to update samples dict given coordinates and a vcf
def update_sample_evidence_dict(sample_evidence_dict, vcf_reader, chrom, pos, ref, alt, study_ids, variant_types):
    '''update sample evidence for a given allele'''

    found_corr_record = False

    if sum([1 for _ in vcf_reader.fetch(chrom, pos - 1, pos)]) != 1:
        print(f"\t[WARNING] more than one record found at pos {pos}")

    for _, record in enumerate(vcf_reader.fetch(chrom, pos - 1, pos)):
        print(
            f"\tparsing samples for row {record.chrom}, {record.pos}, {record.ref}, {record.alts}"
        )

        assert not found_corr_record, \
            "more than 1 record corresponds to the same ref, alt"

        assert ref == record.ref, f"expected ref {record.ref}, got {ref}"
        if alt not in record.alts:
            print(
                f"\tSkipping {chrom}-{pos}-{ref}-{alt}, expected one of {record.alts} alts, got {alt}"
            )
            continue

        # populate patient/sample dictionary
        found_corr_record = True
        id_count = 0

        for sample, genotype in record.samples.items():
            gt_values = genotype["GT"]

            # Check if any allele is non-zero
            if any(int(allele) > 0 for allele in gt_values):
                id_count += 1
                sample_evidence_dict[sample] = {
                    "study_ids": list(study_ids),
                    "variant_types": list(variant_types)
                }

        return id_count


# get focus allele count for cohort allele frequency given a vcf and coordinates
def get_focus_allele_count(vcf_reader, chrom, pos, ref, alt):
    '''calculates focus allele count for a given allele'''

    for _, record in enumerate(vcf_reader.fetch(chrom, pos - 1, pos)):
        print(f"\tCAF: getting allele counts for {record.chrom}, " + \
              f"{record.pos}, {record.ref}, {record.alts}")

        # check correct row
        assert ref == record.ref, f"expected ref {record.ref}, got {ref}"
        if alt not in record.alts:
            print(f"\tSkipping {chrom}-{pos}-{ref}-{alt}, expected one of " + \
                f"{record.alts} alts, got {alt}")
            continue

        # get number of alleles across all patients
        focus_allele_count = sum(
            [sum(genotype["GT"]) for _, genotype in record.samples.items()]
        )

    return focus_allele_count


# create a structured caf dict
def create_caf_dict(
    allele_id,
    gnomad_expression,
    focus_allele,
    focus_allele_count,
    locus_allele_count,
    ancillary_results=None,
):
    '''create cohort allele frequency dict'''
    allele_frequency = focus_allele_count * 1.0 / locus_allele_count

    caf_dict = {
        "id": allele_id,
        "type": "CohortAlleleFrequency",
        "label": f"Overall Cohort Allele Frequency for {allele_id} ({gnomad_expression})",
        "derivedFrom": {
            "id": f"AnVIL_1000G_PRIMED-data-model",
            "type": "DataSet",
            "label": f"AnVIL 1000G PRIMED Data Model",
            "version": "Last Updated 3/12/2024",
        },
        "focusAllele": focus_allele,
        "focusAlleleCount": focus_allele_count,
        "locusAlleleCount": locus_allele_count,
        "alleleFrequency": allele_frequency,
        "cohort": {
            "id": "all_3202_samples",
            "label": "High-coverage sequencing of 3202 samples",
        },
    }

    if ancillary_results is not None:
        caf_dict["ancillaryResults"] = ancillary_results

    return caf_dict

# data stores
sample_evidence_dict = defaultdict(list)
variant_types_set = set()
caf_dicts = []
metakb_api_hits = []

# parse variant matches to save metakb study data
for file_path, matches in matches_per_file.items():
    vcf_reader = pysam.VariantFile(file_path)

    print(truncate(file_path, 0, 47))

    # for each variant per file
    for allele_id, allele_info in matches.items():
        sample_dict = defaultdict(list)

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
        else:
            metakb_api_hits.append(allele_id)
            study_ids = metakb_response["study_ids"]

            assert study_ids == [
                study["id"] for study in metakb_response["studies"]
            ], "study ids in wrong order to add variant types"
            variant_types = [
                study["qualifiers"]["alleleOrigin"]
                for study in metakb_response["studies"]
            ]
            variant_types_set.update(variant_types)

            for study in metakb_response["studies"]:
                print(f"\t{study['type']} ({study['id']}): {study['description']}")
            print()

        id_count = update_sample_evidence_dict(sample_dict, \
                  vcf_reader, chrom, pos, ref, alt, study_ids, variant_types)

        print(f"\tNumber of matching samples: {id_count}\n")

        # get cohort allele counts
        focus_allele_count = get_focus_allele_count(vcf_reader, \
                           chrom, pos, ref, alt)
        print("\tFocus allele count: ", focus_allele_count)

        locus_allele_count = num_samples * 2
        allele_freq = focus_allele_count * 1.0 / (num_samples * 2)
        focus_allele = allele_info[VRS_OBJECT] if VRS_OBJECT in allele_info else allele_id

        # create ancillary_results for caf
        ancillary_results = {
            "patient_matches": id_count,
            "sample_dict": sample_dict
        }

        if metakb_response is not None:
            ancillary_results["metakb_dict"] = metakb_response

        # create and save caf dict
        caf_dict = create_caf_dict(
            allele_id,
            gnomad_expr,
            focus_allele,
            focus_allele_count,
            locus_allele_count,
            ancillary_results,
        )
        print(f"\tAdding {caf_dict['label']}\n")
        caf_dicts.append(caf_dict)

    print("~~~~~~~~~~~~~\n")

# output stats on cache to metakb api hits
metakb_cache_hits = sum(len(e) for e in matches_per_file.values())
print(f"ratio of metakb cache hits to api hits: {len(metakb_api_hits)}/{metakb_cache_hits}")
print("Note: some metakb cache hits are not api hits because some VRS IDs may have matched molecular profiles", \
     "that have no corresponding study ID in metakb because it did not pass metakb's normalizers")

# write Cohort Allele Frequency (CAF) objects to file, then read from it
caf_dir = os.path.expanduser(f"{BASE_DIR}/state")
os.makedirs(caf_dir, exist_ok=True)
file_name = f'caf_objects_{TIMESTAMP}.json'

with open(f'{caf_dir}/{file_name}', 'w') as file:
    json.dump(caf_dicts, file)

print(f"\nwrote cohort allele frequency objects to {caf_dir}/{file_name}!")
