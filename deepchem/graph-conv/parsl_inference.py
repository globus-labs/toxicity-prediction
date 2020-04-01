from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider
from parsl.config import Config
from parsl import python_app
import parsl

from concurrent.futures import as_completed
from multiprocessing import Pool
from datetime import datetime
from rdkit import Chem
import pandas as pd
import numpy as np
import argparse
import logging
import json
import os


def is_smiles_valid(smiles: str) -> bool:
    """Check whether a SMILES string is valid"""
    return Chem.MolFromSmiles(smiles) is None

@python_app()
def inference_function(smiles, **other_cols):
    """Run inference on a list of smiles
    
    Uses multi-processing for intra-node parallelism"""
    # Launch the process pool if this is the first invocation
    #  Note: The pool will stay alive until the host process dies
    #   OK for HPC (host dies when job completes) but be very careful
    #   running this function on persistant servers.
    #global pool 
    import os
    core_count = len(os.sched_getaffinity(0))
    # I use the affinity rather than `os.cpu_count()` to work with aprun's
    #  protocol for specifying the affinity of each MPI PE and all its 
    #  child processes (including those spawned by multiprocessing)
    if 'pool' not in globals():
        from multiprocessing import Pool
        pool = Pool(core_count)
    
    # Measure the start time and record host name
    from datetime import datetime
    from platform import node
    start_time = datetime.utcnow().isoformat()
    hostname = node()
    
    # Pull in the inference function and run it
    from gctox.model import invoke_model
    from gctox.features import compute_features
    import numpy as np
    n_splits = min(core_count * 2, len(smiles))
    chunks = np.array_split(smiles, n_splits)
    feats = np.concatenate(pool.map(compute_features, chunks))
    result = invoke_model(feats, smiles)
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


def write_results(result: dict, path): 
    """Write result of inference function to disk
    
    Appends to existing data file
    """
    # Parse the data
    data = pd.DataFrame(result['result'])
    exists = os.path.isfile(path)

    # Get the runtime and save it
    data['runtime'] = (datetime.fromisoformat(result['end']) - datetime.fromisoformat(result['start'])).total_seconds()
    data['start_time'] = result['start']
    data['end_time'] = result['end']
    data['host'] = result['hostname']
    data['core_count'] = result['core_count']

    # Save the result to disk
    data.to_csv(path, mode='a', header=not exists, index=False)

        
if __name__ == "__main__":
    
    # Get the path to the file to screen
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", help="Number of SMILES strings to send per batch",
                       default=1024, type=int)
    parser.add_argument("--local-workers", help="Number of workers to use for SMILES validation",
                       default=None, type=int)
    parser.add_argument("file", help="Path to the input file")
    args = parser.parse_args()
    
    # Set up the logging
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
    
    # Reading the input data file
    database = pd.read_csv(args.file, header=None, delim_whitespace=True)
    logging.info(f'Loaded {len(database)} molecules from {args.file}')
    
    # Sanitize the input data
    database.rename(columns={0: 'smiles', 1: 'identifier'}, inplace=True)
    with Pool(args.local_workers) as p:
        database['invalid'] = p.map(is_smiles_valid, database['smiles'])
    database.query('not invalid', inplace=True)
    logging.info(f'Found {len(database)} valid SMILES')
    
    # Define the Parsl configuration
    config = Config(
        executors=[
            HighThroughputExecutor(
                address="localhost",
                label="htex",
                max_workers=1,
                provider=LocalProvider(
                    init_blocks=1,
                    max_blocks=1
                ),
            ),
        ],
        strategy=None,
    )
    parsl.load(config)
    logging.info('Configured Parsl')
    
    # Run the inference function
    futures = [inference_function(chunk['smiles'].tolist(), identifier=chunk['identifier'].tolist())
               for chunk in np.array_split(database, len(database) // args.batch_size)]
    logging.info(f'Submitted {len(futures)} tasks')
    
    # Prepare the output file
    output_path = f'{os.path.basename(args.file).split(".")[0]}.csv'
    if os.path.isfile(output_path):
        logging.info(f'Removing old file: {output_path}')
        os.unlink(output_path)
    
    # Receiving results
    for future in as_completed(futures):
        write_results(future.result(), output_path)
        
    logging.info(f'Finished writing to {output_path}')
