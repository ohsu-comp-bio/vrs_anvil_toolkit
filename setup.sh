#!/bin/bash

# download seqrepo data locally
SEQREPO_ROOT=~/seqrepo
mkdir $SEQREPO_ROOT

# for local use
seqrepo --root-directory $SEQREPO_ROOT pull -i 2016-08-28
# # TODO: to replace on Terra / with enough space
# seqrepo --root-directory $SEQREPO_ROOT pull --update-latest

# store new directory
echo "SEQREPO_LATEST_DIR=$(
    find $SEQREPO_ROOT -depth 1 -type d -exec stat -f '%m %N' {} + \
    | sort -rn \
    | head -n 1 \
    | cut -f2- -d' '
)" >> .env