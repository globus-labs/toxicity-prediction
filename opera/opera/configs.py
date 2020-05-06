from parsl import HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.config import Config
from parsl.launchers import AprunLauncher
from parsl.providers import LocalProvider, CobaltProvider
import os

# Read in the environmental variables
with open(os.path.join(os.path.dirname(__file__), '..', 'set_envs.sh')) as fp:
    envs = fp.read()

local_config = Config(
    executors=[
        HighThroughputExecutor(
            address="localhost",
            label="htex",
            max_workers=1,
            prefetch_capacity=1,
            provider=LocalProvider(
                init_blocks=1,
                max_blocks=1
            ),
        ),
    ],
    strategy=None,
)

config = Config(
        retries=1,
        usage_tracking=True,
    executors=[
        HighThroughputExecutor(
            address=address_by_hostname(),
            label="htex",
            max_workers=4,
            prefetch_capacity=1,
            provider=CobaltProvider(
                queue='CVD_Research',
                account='CVD_Research',
                launcher=AprunLauncher(overrides="-d 256 --cc depth -j 4"),
                walltime='3:00:00',
                nodes_per_block=4,
                init_blocks=1,
                min_blocks=1,
                max_blocks=4,
                scheduler_options='#COBALT --attrs enable_ssh=1',
                worker_init=f'''
module load miniconda-3
module load java
source activate /home/lward/exalearn/covid/toxicity-prediction/opera/env
export KMP_AFFINITY=disabled
which python

# Set the environment variables
{envs}
''',
                cmd_timeout=120,
            ),
        ),
    ],
    strategy=None,
)

