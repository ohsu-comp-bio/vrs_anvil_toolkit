import ast
import json
import logging
import multiprocessing
import os
import pickle
import pysam
import requests

from biocommons.seqrepo import SeqRepo
from datetime import datetime
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator

def seqrepo_dir():
   with open(".env") as f:
      for line in f:
         # Ignore comments and empty lines
         if line.strip() and not line.strip().startswith('#'):
            key, value = line.strip().split('=', 1)
            if key == "SEQREPO_ROOT":
               return value + "/latest"

def annotate_vcf(path, require_validation=False):
    '''param stem: path of input vcf file'''
    stem = path.replace(".vcf", "")
    
    input_vcf = path
    # output_vcf = f"{stem}.output.vcf.gz"
    output_vcf = ""
    output_pkl = f"{stem}-vrs-objects.pkl"

    vcf_annotator = VCFAnnotator(seqrepo_root_dir=seqrepo_dir())
    vcf_annotator.tlr.rle_seq_limit = None
    vcf_annotator.annotate(vcf_in=input_vcf, vcf_out=output_vcf, \
        vrs_pickle_out=output_pkl, require_validation=require_validation)
    
    return output_vcf, output_pkl

def print_dict(d, indent=2):
    """pretty print object as json"""
    print(json.dumps(d, indent=indent))

def print_percent(a, b):
    "pretty print percentages"
    print(f"{a}/{b} = {(100.0*a/b):.1f}%")

def annotate_vcf(input_vcf, output_vcf, output_pkl, seqrepo_root_dir, require_validation=True, rle_seq_limit=50):
    '''param stem: path of input vcf file'''
    vcf_annotator = VCFAnnotator(seqrepo_root_dir=seqrepo_root_dir)
    vcf_annotator.tlr.rle_seq_limit = rle_seq_limit
    vcf_annotator.annotate(vcf_in=input_vcf, vcf_out=output_vcf, \
        vrs_pickle_out=output_pkl, require_validation=require_validation)
    
    return output_vcf, output_pkl

def unpickle(file_name):
    """Unpickle vrs objects to single dict"""
    with open(file_name, 'rb') as f:
        vrs_objects = pickle.load(f)
        for k, v in vrs_objects.items():
            yield k, v

def meta_kb(id, recent=True, log=False):
    """Query metakb using vrs object"""
    # k, allele_dict = item
    
    # if translator is not None:
    #     if log: print(f"by {fmt}...")
    #     allele = models.Allele(**allele_dict)
    #     _id = translator.translate_to(allele, fmt)
    # else:
    #     if log: print("by vrs id...")
    #     _id = allele_dict["id"]
        
        
    if recent:
        if log: print("recent elasticbeanstalk api (VRS 2.0 models)")
        response = requests.get("http://metakb-dev-eb.us-east-2.elasticbeanstalk.com" \
                                f"/api/v2/search/studies?variation={id}")
    else:
        if log: print("old api (VRS 1.3 models)")
        response = requests.get("https://dev-search.cancervariants.org" \
                                f"/api/v2/search?variation={id}&detail=false")
    
    response_json = response.json()
    
    if response_json['warnings'] == []:
        return (id, response_json)
    else:
        if log: print(response_json['warnings'])

def num_variants(input_vcf):
    # get total num_variants
    vcf_reader = pysam.VariantFile(open(input_vcf, 'r'))
    return sum(1 for _ in vcf_reader)

def allele_dict_to_hgvs(allele_dict):
    data_proxy = SeqRepoDataProxy(SeqRepo("/Users/quinnwaiwong/seqrepo/latest"))
    translator = AlleleTranslator(data_proxy)
    translator.rle_seq_limit=None
    fmt = "hgvs"
    allele = models.Allele(**allele_dict)
    return translator.translate_to(allele, fmt)

def parallelize(vrs_decorator, vrs_objects, worker_count, progress_interval=500, limit=None):
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
    
    return results