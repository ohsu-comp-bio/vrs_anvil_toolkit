# Cohort Allele Frequency Generation

### Description
Create a cohort allele frequency object for a given variant, subsettable by participant list, GREGOR-formatted phenotypes, etc.

### General Prerequisites
- Variant ID of interest
- VCF Path to file
- Access to phenotypes table either through Terra (default) or as a local file (structured according to the [GREGOR data model](https://gregorconsortium.org/data-model))

### Use Cases
1. Given a variant ID and VCF path, get the allele frequency for the entire cohort
   - Get VCF row corresponding to variant ID using a variant -> VCF row index
   - Get phenotypes corresponding to each participants using the phenotypes by patient table
   - Aggregate counts for participants using their genotypes
   - Create CAF object using counts

2. Given a variant ID, VCF path, **and participant list**, get the allele frequency for a subset of participants (subcohort)
   - Same as 1, just subsetted on a participant list

3. Given a variant ID, VCF path, **and phenotype**, get the allele frequency for the cohort conditional on the phenotype
   - Same as 1, but only increase the counts for the variant of interest if a given patient has the specified phenotype

### Arguments
 - `variant_id` (String): variant ID of interest (VRS ID)
 - `vcf_path` (String): path to VCF file
 - `phenotype_table` (String, optional): where to pull phenotype information from. Defaults to None.
 - `participant_list` (List of Strings, optional): Subset of participants to use. Defaults to None.
 - `phenotype` (String, optional): Specific phenotype to subset on. Defaults to None.

### Caveats
- For multiple alleles, the cohort allele frequency returned is based only on the position and not on the state. In other words, all alleles are on a given variant are handled together.
- For chromosomes with ploidy of 1 (mitochondrial calling or sex chromosomes), focus allele counts (AC) and locus allele counts (AN) can have a maximum value of 1. Focus allele counts are 1 when the genotype has at least a single allele match (0/1, 1/1, or 1) otherwise it is none.
