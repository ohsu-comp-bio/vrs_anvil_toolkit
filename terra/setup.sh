#!/bin/bash

pip install .

# setup reference sequence database (seqrepo)
SEQREPO_ROOT="$HOME/seqrepo"
SEQREPO_TARBALL="seqrepo.tar.gz"

gsutil cp "$WORKSPACE_BUCKET/$SEQREPO_TARBALL" $HOME
echo "creating seqrepo from tar..."
tar -xzf "$HOME/$SEQREPO_TARBALL" --directory=$HOME
seqrepo --root-directory $SEQREPO_ROOT update-latest
echo "removing tarball..."
rm "$HOME/$SEQREPO_TARBALL"
echo "done"
