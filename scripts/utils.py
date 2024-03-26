import ast
import json
import logging
import multiprocessing
import os
import pickle
import subprocess
import pysam
import requests

from biocommons.seqrepo import SeqRepo
from datetime import datetime
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator
from pathlib import Path
from vrs_anvil import seqrepo_dir

# TODO: remove after adding MetaKB API query functionality


# get vrs ids
def translate(gnomad_expr):
    data_proxy = SeqRepoDataProxy(SeqRepo(seqrepo_dir()))
    translator = AlleleTranslator(data_proxy)
    allele = translator._from_gnomad(gnomad_expr)
    return (gnomad_expr, dict(allele))


def annotate_vcf(path, require_validation=False):
    """param stem: path of input vcf file"""
    stem = path.replace(".vcf", "")

    input_vcf = path
    # output_vcf = f"{stem}.output.vcf.gz"
    output_vcf = ""
    output_pkl = f"{stem}-vrs-objects.pkl"

    vcf_annotator = VCFAnnotator(seqrepo_root_dir=seqrepo_dir())
    vcf_annotator.tlr.rle_seq_limit = None
    vcf_annotator.annotate(
        vcf_in=input_vcf,
        vcf_out=output_vcf,
        vrs_pickle_out=output_pkl,
        require_validation=require_validation,
    )

    return output_vcf, output_pkl


def print_dict(d, indent=2):
    """pretty print object as json"""
    print(json.dumps(d, indent=indent))


def print_percent(a, b):
    "pretty print percentages"
    print(f"{a}/{b} = {(100.0*a/b):.1f}%")


def annotate_vcf(
    input_vcf,
    output_vcf,
    output_pkl,
    seqrepo_root_dir,
    require_validation=True,
    rle_seq_limit=50,
):
    """param stem: path of input vcf file"""
    vcf_annotator = VCFAnnotator(seqrepo_root_dir=seqrepo_root_dir)
    vcf_annotator.tlr.rle_seq_limit = rle_seq_limit
    vcf_annotator.annotate(
        vcf_in=input_vcf,
        vcf_out=output_vcf,
        vrs_pickle_out=output_pkl,
        require_validation=require_validation,
    )

    return output_vcf, output_pkl


def unpickle(file_name):
    """Unpickle vrs objects to single dict"""
    with open(file_name, "rb") as f:
        vrs_objects = pickle.load(f)
        for k, v in vrs_objects.items():
            yield k, v


def get_num_variants(input_vcf):
    # get total num_variants
    return int(
        subprocess.run(
            f"grep -v '^#' {input_vcf} | wc -l",
            stdout=subprocess.PIPE,
            shell=True,
            check=True,
            text=True,
        ).stdout
    )


def parallelize(
    vrs_decorator, vrs_objects, worker_count, progress_interval=500, limit=None
):
    """harvest data from service"""

    manager = multiprocessing.Manager()
    results = manager.list()

    with multiprocessing.Pool(worker_count) as pool:
        # call the function for each item in parallel
        c = 0
        print(datetime.now().isoformat(), c)

        for result in pool.imap(vrs_decorator, vrs_objects):
            c += 1
            if result:
                results.append(result)
            if c == limit:
                break
            elif c % progress_interval == 0:
                print(datetime.now().isoformat(), c)

    return list(results)


def truncate(s, first_few, last_few):
    "truncate string printing only first_few and last_few characters"
    return f"{s[:first_few]}...{s[-last_few:]}"
