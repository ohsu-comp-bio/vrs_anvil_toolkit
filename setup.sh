#!/bin/bash

# setup conda environment for use
conda env create -f env.yaml

# download seqrepo data locally
export SEQREPO_ROOT=~/seqrepo
echo $SEQREPO_ROOT
mkdir $SEQREPO_ROOT
# seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

# for local use
seqrepo --root-directory $SEQREPO_ROOT pull -i 2016-08-28