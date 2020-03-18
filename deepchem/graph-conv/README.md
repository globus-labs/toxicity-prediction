# Graph Convolution

The "graph convolution" network from Deepchem was reported to be the best performer for Tox21. 
This directory contains code needed to train, save and run inference at scale with this model.

## Installation

First run the model training code, which should take a few minutes on a laptop with: `python train_model.py`.

Then, install the utilities for running inference on the model with `pip install -e .`.

## Running Inference

The code is designed to provide a single function to run inference on model that we can plug into many workflow 
engines. An example for running from the command line is in `python run_inference.py`

### FuncX

I have registered a function in FuncX for the inference. It must be run on an endpoint configured to use a 
Python environment with the packages described above. See `funcx_inference.py` for an example on how to run 
the predictions. 

Re-register the function by running `register_model.py`.
The function we register to FuncX is designed for high-performance inference on a many-core processor. 
We optimize performance by performing feature generation in parallel with Python's multiprocessing
and using Tensorflow's innate parallelism for the model inference (as opposed to running many models in parallel).
Python startup times are minimized by separating the feature generation from model inference in the 
module so that only one process per node loads the model from disk.
