import parsl

from concurrent.futures import as_completed
from datetime import datetime
import pandas as pd
import numpy as np
import argparse
import logging
import os

from opera.configs import local_config as config
from opera.inference import inference_function


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
    parser.add_argument("file", help="Path to the input file")
    args = parser.parse_args()
    
    # Set up the logging
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
    
    # Reading the input data file
    database = pd.read_csv(args.file, header=None, delim_whitespace=True)
    logging.info(f'Loaded {len(database)} molecules from {args.file}')
    
    # Sanitize the input data
    database.rename(columns={0: 'smiles', 1: 'identifier'}, inplace=True)
    
    # Define the Parsl configuration
    parsl.load(config)
    logging.info('Configured Parsl')
    
    # Run the inference function
    futures = [inference_function(chunk['smiles'].tolist(), identifier=chunk['identifier'].tolist())
               for chunk in np.array_split(database, len(database) // args.batch_size)]
    logging.info(f'Submitted {len(futures)} tasks')
    
    # Prepare the output file
    output_path = os.path.join('screening-results', f'{os.path.basename(args.file).split(".")[0]}.csv')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if os.path.isfile(output_path):
        logging.info(f'Removing old file: {output_path}')
        os.unlink(output_path)
    
    # Receiving results
    for future in as_completed(futures):
        write_results(future.result(), output_path)
        
    logging.info(f'Finished writing to {output_path}')
