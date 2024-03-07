'''
Two goals
- given hgvs IDs, ensure ~100% have hits in metakb
- convert hgvs -> vrs IDs to try to match % hits in metakb as above
'''

import ast
from datetime import datetime
import logging

import pickle
import pysam
import os
import sys
from tqdm import tqdm
import multiprocessing

from biocommons.seqrepo import SeqRepo
from functools import partial
from ga4gh.core._internal.identifiers import ga4gh_identify
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator, Translator
from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator

from utils import seqrepo_dir, meta_kb, parallelize, print_dict, print_percent

if __name__ == '__main__':
    # load in hgvs and convert to vrs ids
    with open("data/genomic_labels.txt", "r") as file:
        hgvs_ids = [line.strip() for line in file.readlines()]

    
    # setup translator

    data_proxy = SeqRepoDataProxy(SeqRepo(seqrepo_dir())) # TODO: not working atm
    
    # seqrepo_rest_service_url = "https://services.genomicmedlab.org/seqrepo"
    # data_proxy = SeqRepoRESTDataProxy(base_url=seqrepo_rest_service_url)

    translator = AlleleTranslator(data_proxy)
    fmt = 'hgvs'
    # translator.rle_seq_limit=None
    # test_hgvs_id = "NC_000013.11:g.32936732C>G"
    # allele1 = translator.translate_from(test_hgvs_id,'hgvs')
    # print(ga4gh_identify(allele1))
    # vrs_objs = [translator._from_hgvs(hgvs_id) for hgvs_id in hgvs_ids]

    alleles = [translator.translate_from(hgvs_id, fmt) for hgvs_id in hgvs_ids[:10]]
    vrs_ids = [ga4gh_identify(allele) for allele in alleles]
    print(list(zip(hgvs_ids, vrs_ids)))
    
    
    hits = []
    hgvs_evidences = []
    hgvs_misses = []
    for label in hgvs_ids:
        evidence = meta_kb(label)
        if evidence is not None:
            hgvs_evidences.append(evidence)
        else:
            hgvs_misses.append(label)

    vrs_hits = parallelize(meta_kb, vrs_ids)

    print(f"Hit percentage:")
    print_percent(len(hgvs_evidences), len(hgvs_ids))
    print(f"Misses: {hgvs_misses}")
    print("metakb: VRS ID hits/HGVS hits")
    print_percent(len(vrs_hits), len(hgvs_evidences))


    # for vrs_id in tqdm(vrs_ids): 
        # potential_hit = meta_kb(vrs_id)
        # print(vrs_id, potential_hit)
        # if potential_hit is not None:
        #     hits.append(potential_hit)


    exit()
            
    for vrs_dict in vrs_dicts:
        hits = []

        # parallelize(meta_kb_by_hgvs, vrs_objects=vrs_dict)
        
        for vrs_obj in tqdm(vrs_dict.items()): 
            try:
                potential_hit = meta_kb(vrs_obj)
            except:
                continue
            if potential_hit:
                print(f"\n ~~~~~~~~ hit! {potential_hit} ~~~~~~~~~~ \n")
                hits.append(potential_hit)