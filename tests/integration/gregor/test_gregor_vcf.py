from collections import defaultdict


def test_read_vcf_associate_participants():
    """Read a VCF, and associate participants with values, decipher values."""
    vcf_fields = "CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT".split()

    # workspace: https://anvil.terra.bio/#workspaces/gregor-dcc/GREGOR_COMBINED_CONSORTIUM_U03/data
    # downloaded from https://storage.googleapis.com/fc-secure-0392c60d-b346-4e60-b5f1-602a9bbfb0e1/gregor_joint_callset_022824/chr_vcfs/gregor_consortium_dataset_u03_v2.chrY.vcf.gz
    samples_with_values = defaultdict(int)
    values = defaultdict(int)
    with open('gregor_consortium_dataset_u03_v2.chrY.vcf', 'r') as file:
        for line in file.readlines():
            # Process each line
            if line.startswith('#CHROM'):
                line = line.replace('#CHROM', 'CHROM')
                header = line.strip().split('\t')
            elif line.startswith('#'):
                continue
            else:
                # create a zip of the header and the line
                line = dict(zip(header, line.strip().split('\t')))
                for key, value in line.items():
                    # skip the standard VCF fields
                    if key in vcf_fields:
                        continue
                    # skip the values that are not populated ?
                    if value == './.:.:.:.':
                        continue
                    # TODO - where are documentation for values?
                    # e.g.: GT:AD:GQ:RGQ, GT:AD:FT:GQ:RGQ
                    # `0/0:.:40:.` vs `1/1:0,16:45:123`
                    # from header

                    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                    # see https://samtools.github.io/hts-specs/VCFv4.1.pdf

                    # In Variant Call Format (VCF) files, GT stands for genotype, which is encoded as allele values separated by a slash (/) or vertical pipe (|):
                    # 0: The reference base
                    # 1: The first entry in the ALT column
                    # 2: The second allele listed in ALT
                    # Forward slash (/): Indicates that no phasing information is available
                    # Vertical pipe (|): Indicates that the genotype is phased

                    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">

                    # see https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format#:~:text=Allele%20depth%20(AD)%20and%20depth,each%20of%20the%20reported%20alleles.
                    # In a Variant Call Format (VCF) file, AD stands for Allele Depth, which is the number of reads that support each allele. The AD field is an array that includes the reference allele as the first entry, and the remaining entries are for each alternate allele at that locus

                    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">

                    # see https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format#:~:text=PL%20is%20calculated.-,GQ,how%20they%20should%20be%20used.
                    # In Variant Call Format (VCF), GQ stands for Genotype Quality and is a measure of how confident a genotype assignment is. GQ is calculated by taking the difference between the PL of the second most likely genotype and the PL of the most likely genotype. The PLs are normalized so that the most likely PL is always 0, so the GQ is usually equal to the second smallest PL. However, the GQ is capped at 99, so if the second most likely PL is greater than 99, the GQ will be 99.

                    ##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">

                    # see https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized
                    # Reference Genotype Quality (RGQ) -- The phred-scaled confidence that the reference genotypes are correct.  A higher score indicates a higher confidence.  For more information on RGQ, please see the GQ documentation, but note that RGQ applies to the reference, not the variant.  For more information on interpreting phred-scaled values, please see Phred-scaled quality scores.

                    # update a counter
                    samples_with_values[key] += 1
                    values[value] += 1

    max_sample = max(samples_with_values, key=samples_with_values.get)
    max_sample_count = samples_with_values[max_sample]

    # Print the sample and its hit count
    print(f"Sample with the greatest hit: {max_sample} with {max_sample_count} hits")

    max_value = max(values, key=values.get)
    max_value_count = values[max_value]

    print(f"Value with the greatest hit: {max_value} with {max_value_count} hits")

    # >>> Sample with the greatest hit: GSS193307 with 272090 hits
    assert False, "TODO - where are documentation for values?"
