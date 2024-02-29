from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from utils import *

from datetime import datetime
import os
import subprocess

def calculate_gnomad_expressions(input_vcf, alt=True):
    if alt: command = f"bcftools query -f '%CHROM-%POS-%REF-%ALT\n' {input_vcf}"
    else: command = f"bcftools query -f '%CHROM-%POS-%REF-%REF\n' {input_vcf}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    # Iterate over the output of bcftools and yield each gnomAD expression
    for line in process.stdout:
        yield line.strip()



if __name__ == "__main__":
    input_vcf = "vcf/NA12878.chr1.100.000.vcf"
    # input_vcf = "vcf/1kGP.chr1.1000.vcf"
    stem = input_vcf.replace('.vcf', '')
    output_pkl = f"{stem}.impl.2-vrs-objects.pkl"

    gnomad_generator_alt = calculate_gnomad_expressions(input_vcf)
    gnomad_generator_ref = calculate_gnomad_expressions(input_vcf, alt=False)


    gnomad_generators = [gnomad_generator_ref, gnomad_generator_alt]
    manager = multiprocessing.Manager()
    allele_dicts_list = [manager.dict() for _ in range(len(gnomad_generators))]
    worker_count = os.cpu_count()
    progress_interval = 10_000

    for allele_dict, gnomad_generator in zip(allele_dicts_list, gnomad_generators):
        manager = multiprocessing.Manager()

        with multiprocessing.Pool(worker_count) as pool:
            # call the function for each item in parallel
            c = 0
            print(datetime.now().isoformat(), c)

            for result in pool.imap(translate, gnomad_generator):
                c += 1
                if result:
                    allele_dict[result[0]] = result[1]
                elif c % progress_interval == 0:
                    print(datetime.now().isoformat(), c)

    # save to pkl
    combined_allele_dicts = {}
    for shared_dict in allele_dicts_list:
        combined_allele_dicts.update(shared_dict)

    with open(output_pkl, 'wb') as file:
        pickle.dump(combined_allele_dicts, file)

    # with open(output_pkl, 'rb') as file:
    #     ting = pickle.load(file)
    #     for k,v in ting.items():
    #         print(k,v)
    #         print(type(v))
    #         break
        
    #     print(len(ting))