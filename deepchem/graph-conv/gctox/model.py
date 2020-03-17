"""Function for running inference"""
from .features import compute_features
from typing import List
from rdkit import Chem
import numpy as np
import json
import os

# Load in the model
import tensorflow as tf
np.random.seed(123)

tf.set_random_seed(123)
import deepchem as dc
from deepchem.data.datasets import NumpyDataset
from deepchem.models.tensorgraph.models.graph_models import GraphConvModel

model_dir = os.path.join(os.path.dirname(__file__), "..", "model")

# Create the featurizer and transformer
with open(os.path.join(model_dir, '..', 'tasks.json'), 'r') as fp:
    tasks = json.load(fp)

# Batch size of models
model = GraphConvModel(12, mode='classification', model_dir=model_dir, batch_size=128)
model.restore()

# Make the inference functions
def invoke_model(feats: np.array, smiles: List[str]) -> [dict]:
    """Invoke the model
    
    Args:
        feats (np.array): Features for the model
        smiles ([str]): SMILES 
    Returns:
        ([dict]) Return the data
    """
    # Turn the features into a Numpy dataset
    dataset = NumpyDataset(feats, n_tasks=len(tasks))

    # Run inference
    y_pred = model.predict(dataset)

    # Get the output
    tox_liklihood = y_pred[:, :, 1]
    output = dict(zip(tasks, tox_liklihood.T))
    output['smiles'] = smiles
    return output


def run_inference(smiles: List[str]) -> [dict]:
    """Run inference on the machine learning models

    Args:
        smiles ([str]): List of SMILES to evaluate
    Returns:
        ([dict]) Dictionary of the toxicity liklihoods
    """
    feats = compute_features(smiles)
    return invoke_model(feats, smiles)
