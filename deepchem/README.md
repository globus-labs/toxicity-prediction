# DeepChem-based Models

These directories contain models derived from the [deepchem](http://github.com/deepchem/deepchem/) framework.
We used the examples from the "MoleculeNet" benchmarks performed using deepchem and, specifically, those
with the highest performance on the Tox21 challenge.

The subdirectories each contain code for training each model and a simple function interface for performing
inference using the trained model.

## Installation

We use Anaconda to define the software environment for the deepchem models:

```bash
conda env create --file environment.yml --force
```

This creates an environment named "deepchem" that runs Python 3.7
