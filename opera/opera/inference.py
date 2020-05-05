"""Utilities for running inference on ADMET models"""
from tempfile import TemporaryDirectory
from opera.cdk import is_valid_smiles
from subprocess import Popen, STDOUT
from parsl import python_app
import re
import os

OPERA_PATH = os.environ['OPERA_PATH']
MCR_PATH = os.environ['MCR_PATH']

_atom_regex = re.compile(r'([A-Z][a-z]?|[a-z])')
_cdk_valences = ['H', 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Hg', 'Rn', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al',
                 'Si', 'P', 'S', 'Cl', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Rb', 'Sr', 'In', 'Sn', 'Sb', 'Te', 'I',
                 'Cs', 'Ba', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Fr', 'Ra', 'Cu', 'Mn', 'Co']


def _cdk_ok(smiles: str, cutoff: int = 50) -> bool:
    """Check whether a molecule can be computed with CDK"""

    # Make sure the molecule does not contain deuterium or tritium
    #  I think these mess up some descriptor calculations that have
    #  improper definitions of "heavy atom"
    if "[2H]" in smiles or "[3H]" in smiles:
        return False

    # Check that CDK can parse this SMILES
    if not is_valid_smiles(smiles):
        return False

    atoms = _atom_regex.findall(smiles)
    # Make sure the molecule is not too big
    if len(atoms) > cutoff:
        return False

    # Make sure we know the valences for all elements
    return all(x.capitalize() in _cdk_valences for x in atoms)


@python_app()
def inference_function(smiles, **other_cols):
    """Run inference on a list of smiles using Opera's CLI application"""
    # Measure the start time and record host name
    from datetime import datetime
    from platform import node
    start_time = datetime.utcnow().isoformat()
    hostname = node()
    core_count = len(os.sched_getaffinity(0))

    # Determine which molecules are suitable for the pipeline
    skip = [not _cdk_ok(s) for s in smiles]

    # Make a temporary directory for the results
    with TemporaryDirectory() as td:
        # Store the SMILES string in a result object
        result = {'smiles': smiles}

        # Write the SMILES and an index to disk
        input_path = os.path.join(td, 'input.smi')
        output_path = os.path.join(td, 'output.csv')
        with open(input_path, 'w') as fp:
            for i, (s, f) in enumerate(zip(smiles, skip)):
                if not f:
                    print(f'{s}\t{i}', file=fp)

        # Invoke the opera executable
        command = [os.path.join(OPERA_PATH, 'run_OPERA_P.sh'), MCR_PATH, '-s', input_path,
                   '-o', output_path, '-Tox', '-P', str(core_count)]
        stdout = os.path.join(td, 'stdout')
        with open(stdout, 'w') as so:
            proc = Popen(command, stderr=STDOUT, stdout=so)
            rc = proc.wait()

        # If there is a non-zero error, raise an exception with that content
        if rc != 0:
            with open(stdout) as fp:
                raise ValueError(fp.read())

        # Read in the results from disk
        with open(output_path) as fp:
            # Get the names of columns in the header
            columns = fp.readline().split(",")

            # Prepare the result dictionary to store the result outputs
            for c in columns[1:]:
                result[c] = []

            # Read them into the result dictionary
            for i, f in enumerate(skip):
                if f:
                    for c in columns[1:]:
                        result[c].append(None)
                else:
                    line = fp.readline().strip()
                    data = line.split(",")
                    assert i == int(data[0]), 'Output data is not sorted!'
                    for c, d in zip(columns[1:], data[1:]):
                        result[c].append(d)

        # Add in the extra columns
        result.update(other_cols)
        result['skipped'] = skip

        # Measure the end time
        end_time = datetime.utcnow().isoformat()
        return {
            'start': start_time,
            'result': result,
            'end': end_time,
            'core_count': core_count,
            'hostname': hostname
        }
