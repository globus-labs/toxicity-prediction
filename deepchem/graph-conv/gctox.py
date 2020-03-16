"""Function for running inference"""
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
from deepchem.feat.graph_features import ConvMolFeaturizer
from deepchem.models.tensorgraph.models.graph_models import GraphConvModel

model_dir = os.path.join(os.path.dirname(__file__), "model")

# Create the featurizer and transformer
feat = ConvMolFeaturizer()
with open(os.path.join(model_dir, '..', 'tasks.json'), 'r') as fp:
    tasks = json.load(fp)

# Batch size of models
model = GraphConvModel(12, mode='classification', model_dir=model_dir, batch_size=128)
model.restore()

# Create prediction function
smiles = ['C', 'CC', 'CCC'] 

def run_inference(smiles: List[str]) -> [dict]:
    """Run inference on the machine learning models

    Args:
        smiles ([str]): List of SMILES to evaluate
    Returns:
        ([dict]) Dictionary of the toxicity liklihoods
    """

    # Create the dataset
    mols = [Chem.MolFromSmiles(x) for x in smiles]
    feats = np.array(feat.featurize(mols))
    dataset = NumpyDataset(feats, n_tasks=len(tasks))

    # Run inference
    y_pred = model.predict(dataset)

    # Get the output
    tox_liklihood = y_pred[:, :, 1]
    output = dict(zip(tasks, tox_liklihood.T))
    output['smiles'] = smiles
    return output
