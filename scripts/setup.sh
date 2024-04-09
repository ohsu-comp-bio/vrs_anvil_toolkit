#!/bin/bash

# create venv
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# # download seqrepo data locally
SEQREPO_ROOT=~/seqrepo
echo "SEQREPO_ROOT=$SEQREPO_ROOT" > .env
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

deactivate
