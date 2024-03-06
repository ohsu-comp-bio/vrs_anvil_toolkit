# vrs-python-testing
## Goal
Provide functionality to connect variant-level VCF information to clinical cancer variant evidence. Also specing out the capabilities of variant normalization with `vrs-python` and evidence gathering with `metakb`

## Setup
### Terra
```
cd ~
git clone https://github.com/ohsu-comp-bio/vrs-python-testing.git
cd vrs-python-testing
gsutil cp 1000g_patient_na12878_evidence.ipynb $WORKSPACE_BUCKET/notebooks/
```

Once you open the notebook, it'll be copied down to your persistent disk in `"$HOME/$WORKSPACE/edit`, so you'll need to copy it back to the repo before committing any changes, eg
```
cp `"$HOME/$WORKSPACE/edit/1000g_patient_na12878_evidence.ipynb` $HOME/vrs-python-testing/
```

### Local
```
bash setup.sh
```

## Testing
```
source venv/bin/activate
pytest .
```

## Notes on Usage

Currently defaults to GRCh38 (add an annotate parameter to change it)


### Contributing
This project is open to contributions from the research community. If you are interested in contributing to the project, please contact the project team.
See the [contributing guide](CONTRIBUTING.md) for more information on how to contribute to the project.
