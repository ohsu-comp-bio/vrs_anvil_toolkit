from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from utils import *

import subprocess

def calculate_gnomad_expressions(input_vcf, alt=True):
    if alt: command = f"bcftools query -f '%CHROM-%POS-%REF-%ALT\n' {input_vcf}"
    else: command = f"bcftools query -f '%CHROM-%POS-%REF-%REF\n' {input_vcf}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    # Iterate over the output of bcftools and yield each gnomAD expression
    for line in process.stdout:
        yield line.strip()

input_vcf = "vcf/NA12878.chr1.100.000.vcf"
# input_vcf = "vcf/1kGP.chr1.1000.vcf"
stem = input_vcf.replace('.vcf', '')
output_pkl = f"{stem}.my.impl.no.val-vrs-objects.pkl"

gnomad_generator_alt = calculate_gnomad_expressions(input_vcf)
gnomad_generator_ref = calculate_gnomad_expressions(input_vcf, alt=False)

# get vrs ids
data_proxy = SeqRepoDataProxy(SeqRepo(seqrepo_dir())) # TODO: not working atm
translator = AlleleTranslator(data_proxy)

# allele_dicts = {expr: dict(translator._from_gnomad(expr)) for expr in gnomad_generator}
allele_dicts = {}
for generator in [gnomad_generator_ref, gnomad_generator_alt]:
    for expr in generator:
        result = translator._from_gnomad(expr)
        if result is not None:
            allele_dicts[expr] = dict(result)

# save to pkl
with open(output_pkl, 'wb') as file:
    pickle.dump(allele_dicts, file)

# with open(output_pkl, 'rb') as file:
#     ting = pickle.load(file)
#     for k,v in ting.items():
#         print(k,v)
#         print(type(v))
#         break

# parallelize(vrs_decorator, vrs_objects, worker_count, progress_interval=500, limit=None):
