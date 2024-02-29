#!/bin/bash

# setup python packages and dependencies
pip install -r ~/requirements.txt
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

# setup bcftools
cd ~
curl -LJO https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
BCF_TOOLS_TAR=$(ls -1t bcftools*.tar.bz2 | head -n 1)
tar -xvf $BCF_TOOLS_TAR
rm $BCF_TOOLS_TAR

mv ~/bcftools*/ $BCF_TOOLS_DIR 
cd $BCF_TOOLS_DIR
./configure prefix=$HOME
make
make install
