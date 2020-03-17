from deepchem.feat.graph_features import ConvMolFeaturizer
from typing import List
from rdkit import Chem
import numpy as np

feat = ConvMolFeaturizer()

def compute_features(smiles: List[str]) -> np.array:
    """Compute the features for a list of molecules
    
    Args:
        smiles ([str]): SMILES of molecules to evaluate
    Returns:
        (np.array) Features for each molecule
    """
    # Create the dataset
    mols = [Chem.MolFromSmiles(x) for x in smiles]
    return np.array(feat.featurize(mols))
