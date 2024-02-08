#!/bin/bash

# create python package
pip install -r requirements.txt

# # download seqrepo data locally
SEQREPO_ROOT=~/seqrepo
echo "SEQREPO_ROOT=$SEQREPO_ROOT" > .env
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

cd ~
curl -LJO https://github.com/vcftools/vcftools/tarball/master

VCF_TOOLS_TAR=$(ls -1t vcftools*.tar.gz | head -n 1)
tar -xzvf $VCF_TOOLS_TAR
rm $VCF_TOOLS_TAR

VCF_TOOLS_DIR=$(ls -1td vcftools* | head -n 1)
cd $VCF_TOOLS_DIR
echo "VCF_TOOLS_DIR=$VCF_TOOLS_DIR" >> .env

./autogen.sh
./configure prefix=$HOME
make
make install

export PERL5LIB="$VCF_TOOLS_DIR/src/perl/"
export VCFTOOLS="$VCFTOOLS_DIR/src/cpp/vcftools"
echo "VCF_TOOLS=$VCF_TOOLS" >> .env