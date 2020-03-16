# Graph Convolution

The "graph convolution" network from Deepchem was reported to be the best performer for Tox21. 
This directory contains code needed to train, save and run inference at scale with this model.

## Installation

First run the model training code, which should take a few minutes on a laptop with: `python train_model.py`.

Then, install the utilities for running inference on the model with `pip install -e .`.

## Running Inference

The code is designed to provide a single function to run inference on model that we can plug into many workflow 
engines. An example for running from the command line is in `python run_inference.py`

