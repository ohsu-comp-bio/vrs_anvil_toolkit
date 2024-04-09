#!/bin/bash

# setup python packages and dependencies
pip install -r requirements.txt

SEQREPO_ROOT=~/seqrepo
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest
