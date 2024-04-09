#!/bin/bash

# setup python packages and dependencies
pip install -r requirements.txt

SEQREPO_ROOT=~/seqrepo
BCFTOOLS_DIR=~/bcftools/
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest
