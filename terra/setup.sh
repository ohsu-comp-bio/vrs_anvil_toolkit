#!/bin/bash

# setup python packages and dependencies
pip install -r ~/requirements.txt
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

# setup vcftools
cd ~
curl -LJO https://github.com/vcftools/vcftools/tarball/master

VCF_TOOLS_TAR=$(ls -1t vcftools*.tar.gz | head -n 1)
tar -xzvf $VCF_TOOLS_TAR
rm $VCF_TOOLS_TAR

mv ~/vcftools*/ $VCF_TOOLS_DIR 
cd $VCF_TOOLS_DIR
./autogen.sh
./configure prefix=$HOME
make
make install