#!/bin/bash

# setup python packages and dependencies
pip install -r requirements.txt

SEQREPO_ROOT=~/seqrepo
BCFTOOLS_DIR=~/bcftools/
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

# setup bcftools
if [ bcftools --version &> /dev/null ]; then
    cd ~
    curl -LJO https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
    BCFTOOLS_TAR=$(ls -1t bcftools*.tar.bz2 | head -n 1)
    tar -xvf $BCFTOOLS_TAR
    rm $BCFTOOLS_TAR

    mv ~/bcftools*/ $BCFTOOLS_DIR
    cd $BCFTOOLS_DIR
    ./configure prefix=$HOME
    make
    make install
else
    echo "bcftools already installed, not installing"
fi
