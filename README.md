# Toxicity Model Repository

This repository is a collection of models that predict the toxicity of
molecules from machine learning.
Each subdirectory houses models from different sources.
The subdirectories contain an environment description (preferrably
using an Anaconda YAML file, and further details about where the models
originate and technical descriptions of their outputs.

## Usage

We plan to build a common interface for these models, but 
such efforts have yet to be started.

The interfaces should at least define the ability to render
toxicity estimates from the SMILES string of a molecule.
