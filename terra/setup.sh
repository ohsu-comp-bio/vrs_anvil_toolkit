#!/bin/bash

pip install -r requirements.txt

# setup reference sequence database (seqrepo)
SEQREPO_ROOT="$HOME/seqrepo"
SEQREPO_TARBALL="seqrepo.tar.gz"

mkdir $SEQREPO_ROOT
gsutil cp "$WORKSPACE_BUCKET/$SEQREPO_TARBALL" $HOME
tar -xzf "$HOME/$SEQREPO_TARBALL"
seqrepo --root-directory $SEQREPO_ROOT update-latest
rm "$HOME/$SEQREPO_TARBALL"
