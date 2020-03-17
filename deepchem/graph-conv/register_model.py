from funcx.sdk.client import FuncXClient
import json

def inference_function(smiles):
    """Run inference on a list of smiles
    
    Uses multi-processing for intra-node parallelism"""
    # Set the affinity for this worker.
    #  Needed for Tensorflow, which likes to use all the cores
    global pool 
    import os
    core_count = len(os.sched_getaffinity(0))
    if 'pool' not in globals():
        from multiprocessing import Pool
        pool = Pool(core_count)
    
    # Measure the start time
    from datetime import datetime
    start_time = datetime.utcnow().isoformat()
    
    # Pull in the inference function and run it
    from gctox.model import invoke_model
    from gctox.features import compute_features
    import numpy as np
    n_splits = min(core_count * 2, len(smiles))
    chunks = np.array_split(smiles, n_splits)
    feats = np.concatenate(pool.map(compute_features, chunks))
    result = invoke_model(feats, smiles)
    
    # Measure the end time
    end_time = datetime.utcnow().isoformat()
    return {
        'start': start_time,
        'result': result,
        'end': end_time,
        'core_count': core_count
    }

# Test run
print(inference_function(['C', 'CCCCC']))
    
# Make the client
fxc = FuncXClient()

# Register and save the function
func_uuid = fxc.register_function(inference_function, description="Infer toxicity based on Tox21 with Deepchem's Graph Convolution")
print(f'Registered function as {func_uuid}')
with open('func_uuid.json', 'w') as fp:
    json.dump(func_uuid, fp)
