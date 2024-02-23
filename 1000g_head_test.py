"""
given a vcf file,
- output vrs objects associated with the variants
- convert vrs ids to hgvs expressions
- use hgvs expressions to hit metakb, logging number of hits
"""

import logging
import pickle
import pysam
from time import time

from ga4gh.core._internal.identifiers import ga4gh_identify
from ga4gh.vrs import models
from utils import *

if __name__ == "__main__":
    # input path
    input_vcf_file_name = "1kGP.chr1.1000.vcf"
    input_dir = "vcf"
    
    # output paths
    output_dir = "split"
    os.makedirs(output_dir, exist_ok=True)
    stem = output_dir + "/" + input_vcf_file_name.replace('.vcf', '')
    output_vcf = f"{stem}.vcf.gz"
    output_pkl = f"{stem}-vrs-objects.pkl"
    metakb_output_pkl = f"{stem}-hits.pkl"

    # vrs annotation fields
    seqrepo_root_dir = seqrepo_dir()
    rle_seq_limit = 50 # reference length sequence limit
    require_ref_validation = True # for reference genome

    # metakb settings
    num_ids_limit = None # number of ids to process, None = all
    progress_interval = 100

    # suppress DEBUG logs (1 per variant)
    logger = logging.getLogger("ga4gh.vrs.extras.vcf_annotation")
    logger.setLevel(level=logging.INFO)

    ## store VCF as VRS Alleles
    input_vcf = f"{input_dir}/{input_vcf_file_name}"
        
    print("trying...", input_vcf)
    t = time()
    output_vcf, output_pkl = annotate_vcf(input_vcf, output_vcf, output_pkl, \
        seqrepo_root_dir, require_ref_validation, rle_seq_limit)
    elapsed_time = time()-t
    print(f"annotation: {(elapsed_time):.2f}s")

    # get pickle totals
    allele_dicts = unpickle(output_pkl)

    # get total num_variants
    t = time()
    vcf_reader = pysam.VariantFile(open(input_vcf, 'r'))
    num_variants = sum(1 for _ in vcf_reader)
    elapsed_time = time()-t
    print(f"get num variants: {(elapsed_time):.2f}s")

    print_percent(len(allele_dicts), num_variants)
    print()

    # convert alleles to vrs ids
    vrs_ids = [ga4gh_identify(models.Allele(**allele_dict)) \
            for _, allele_dict in allele_dicts.items()]
            
    # ping metakb
    print("pinging metakb...")
    t = time()
    # number of workers
    worker_count = 4*os.cpu_count()
    hits = parallelize(meta_kb, vrs_ids, worker_count, \
        progress_interval=progress_interval, limit=num_ids_limit)
    print(f"metakb: {(time()-t):.2f} s")

    # save and report hits
    with open(metakb_output_pkl, 'wb') as file:
        pickle.dump(hits, file)
    
    print("\nhits to ids queried...")
    total = num_ids_limit if num_ids_limit else len(vrs_ids)
    print_percent(len(hits), total)

    