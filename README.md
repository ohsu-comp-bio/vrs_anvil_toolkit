<img width="685" alt="image" src="https://github.com/ohsu-comp-bio/vrs-python-testing/assets/47808/909db052-972c-4508-a2f4-8a389de03320">


# VRS AnVIL

## Project Overview

This Python project is designed to process Variant Call Format (VCF) files or other sources of variant information and perform lookup operations on Genomic Variation Representation Service (GA4GH VRS) identifiers. The GA4GH VRS identifiers provide a standardized way to represent genomic variations, making it easier to exchange and share genomic information.

In addition, this project facilitates the retrieval of evidence associated with genomic alleles by leveraging the Genomic Data Representation and Knowledge Base (GA4GH MetaKB) service. GA4GH MetaKB provides a comprehensive knowledge base that links genomic variants to relevant evidence, enabling users to access valuable information about genomic alleles.

## Features

1. **VCF File Processing:**
   - Streamlines reading and parsing of VCF files, to extract relevant genomic information.

2. **GA4GH VRS Identifier Lookup:**
   - Utilizes the GA4GH VRS API to perform lookups for each genomic variation mentioned in the VCF file.
   - Retrieves standardized identifiers for the alleles, enhancing interoperability with GA4GH-compliant systems.
   - GA4GH MetaKB Service Integration:  Utilizes the GA4GH MetaKB service to query and retrieve evidence associated with the specified genomic alleles.
3. **Output Generation:**
   - Generates summary metrics about throughput, errors, and evidence hits and misses
   - Optionally, generates a processed VCF file with additional GA4GH VRS identifiers for each genomic variation.
   - Presents the retrieved evidence in a structured format, including information about studies, publications, and other relevant details.


4. **Error Handling:**
   - Implements robust error handling to address issues like invalid input files, invalid variants, connectivity problems with the GA4GH MetaKB API, and more.

## Getting Started

### Prerequisites

- Python 3.10 or later
- Internet connectivity for setting up dependencies, GA4GH MetaKB lookups, etc.

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/ohsu-comp-bio/vrs-anvil
   cd vrs-anvil
   ```

2. Install dependencies:
   a. for local use
   ```bash
   # install postgresql@14 (required for vrs-python)
   brew install postgresql@14
   bash scripts/setup.sh
   ```
   b. for use on Terra
   ```bash
   SEQREPO_ROOT=~
   bash terra/setup.sh
   ```

### Usage
**Manifest**

The configuration is controlled by a [manifest.yaml](tests/fixtures/manifest.yaml) file. The manifest file specifies the input VCF file(s), the output directory, and other parameters.

**CLI**
```bash
source venv/bin/activate
# navigate to a working directory, with your manifest.yaml file.  Add the VCF urls or file paths to your manifest

# run the vrs_anvil command in the fore ground
vrs_anvil annotate

# run the vrs_anvil command in parallel, one process per VCF file
vrs_anvil annotate --scatter

# run the vrs_anvil command in parallel in the background
nohup vrs_anvil annotate --scatter & # press enter to continue

# get the status of the scatter processes
vrs_anvil ps
```

**Processing VCF Files ([vrs-python](https://github.com/ga4gh/vrs-python))**

vrs-python is a GA4GH GKS package centered around creating Variant Representation specification (VRS) IDs: consistent, globally unique identifiers for variation. Some of its functionality includes variant ID translation and VCF annotation. Used as a dependency in vrs_anvil, it can also be used as a standalone package.

For Python usage, see [vrs_vcf_annotator.py](scripts/vrs_vcf_annotator.py) for an example.

For CLI usage:
```bash
python3 -m ga4gh.vrs.extras.vcf_annotation --vcf_in tests/fixtures/1kGP.chr1.1000.vcf --vcf_out annotated_output.vcf.gz --vrs_pickle_out allele_dicts.pkl --seqrepo_root_dir ~/seqrepo/latest
```

The above is an example using an example vcf. Replace the `--vcf_out` and `vrs_pickle_out` here with your desired output file path, where the output vcf can be BCF (`vcf.gz`) or VCF (`vcf`)

**Terra**
The command line utility supports Google Cloud URIs and running commands in the background to interop with Terra out-of-the-box. This is described in the [CLI usage](#features) above. For an example notebook, see `vrs-anvil-demo.ipynb` on the `vrs-anvil` workspace.

### Contributing

This project is open to contributions from the research community. If you are interested in contributing to the project, please contact the project team.
See the [contributing guide](CONTRIBUTING.md) for more information on how to contribute to the project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.
