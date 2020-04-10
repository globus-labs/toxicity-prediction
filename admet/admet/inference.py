"""Utilities for running inference on ADMET models"""
from parsl import python_app


pool = None
models = None
feat = None


def run_inference(model, chunk):
    return model.predict_proba(chunk)[:, 1]


@python_app()
def inference_function(smiles, model_dir, **other_cols):
    """Run inference on a list of smiles

    Uses multi-processing for intra-node parallelism"""
    # Launch the process pool if this is the first invocation
    #  Note: The pool will stay alive until the host process dies
    #   OK for HPC (host dies when job completes) but be very careful
    #   running this function on persistent servers.
    global pool, models, feat
    import os
    core_count = len(os.sched_getaffinity(0))
    # I use the affinity rather than `os.cpu_count()` to work with aprun's
    #  protocol for specifying the affinity of each MPI PE and all its
    #  child processes (including those spawned by multiprocessing)
    if pool is None:
        from multiprocessing import Pool
        pool = Pool(core_count)

    # Measure the start time and record host name
    from datetime import datetime
    from platform import node
    start_time = datetime.utcnow().isoformat()
    hostname = node()

    # Load models
    from glob import glob
    import pickle as pkl
    if models is None:
        model_files = glob(os.path.join(model_dir, '*.pkl'))
        models = {}
        for path in model_files:
            with open(path, 'rb') as fp:
                models[os.path.basename(path)[:-4]] = pkl.load(fp)

        # Remove the first step from the pipeline (feature generation from the SMILES)
        feat = None
        for m in models.values():
            feat = m.steps.pop(0)[1]

    # Compute features in parallel
    import numpy as np
    n_splits = min(core_count * 4, len(smiles))
    smiles_chunks = np.array_split(smiles, n_splits)
    feature_chunks = pool.map(feat.transform, smiles_chunks)

    # Pull in the inference function and run it
    from functools import partial
    result = {'smiles': smiles}
    for name, model in models.items():
        func = partial(run_inference, model)
        result[name] = np.concatenate(pool.map(func, feature_chunks))
    result.update(other_cols)

    # Measure the end time
    end_time = datetime.utcnow().isoformat()
    return {
        'start': start_time,
        'result': result,
        'end': end_time,
        'core_count': core_count,
        'hostname': hostname
    }
