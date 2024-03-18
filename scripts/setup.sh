#!/bin/bash

# create python package
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# # download seqrepo data locally
SEQREPO_ROOT=~/seqrepo
BCFTOOLS_DIR=~/bcftools/
echo "SEQREPO_ROOT=$SEQREPO_ROOT" > .env
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
    ./configure
    make
    sudo make install
else
    echo "bcftools already installed, not installing"
fi

deactivate
