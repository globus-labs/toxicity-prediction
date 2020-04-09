from typing import Optional, List

import numpy as np
from pybel import readstring
from sklearn.base import BaseEstimator, TransformerMixin


class SubstructureMatching(BaseEstimator, TransformerMixin):
    """Feature vector where each entry indicates the presence of a certain group"""

    def __init__(self, fingerprints: Optional[List[str]] = None):
        """
        Args:
            fingerprints ([str]): List of fingerprints to compute
        """
        if fingerprints is None:
            fingerprints = ['maccs', 'fp4']
        self.fingerprints = fingerprints

        # Compute the size of the fingerprints
        test_smi = readstring('smi', 'C')
        self.fp_size = []
        for f in self.fingerprints:
            fp = test_smi.calcfp(f)
            self.fp_size.append(len(fp.fp) * 32)  # Assumes unsigned ints are

    def fit(self, X, y=None):
        return self

    @property
    def num_features(self):
        return sum(self.fp_size)

    def transform(self, X, y=None):
        # Initialize output array
        output = np.zeros((len(X), sum(self.fp_size)), dtype=np.bool)

        # Compute the fingerprints for each SMILES in the input array
        for i, smi in enumerate(X):
            mol = readstring('smi', smi)

            # Compute each type of fingerprints
            for fi, f in enumerate(self.fingerprints):
                fp = mol.calcfp(f)
                offset = sum(self.fp_size[:fi])
                output[i][np.add(offset, fp.bits).tolist()] = True

        return output
