#!/bin/bash

# create python package
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# # download seqrepo data locally
SEQREPO_ROOT=~/seqrepo
echo "SEQREPO_ROOT=$SEQREPO_ROOT" > .env
mkdir $SEQREPO_ROOT
seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

# get latest directory name
echo "SEQREPO_LATEST_DIR=$(
    find $SEQREPO_ROOT -depth 1 -type d -exec stat -f '%m %N' {} + \
    | sort -rn \
    | head -n 1 \
    | cut -f2- -d' '
)" >> .env

deactivate