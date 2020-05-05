from parsl import HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.config import Config
from parsl.launchers import AprunLauncher
from parsl.providers import LocalProvider, CobaltProvider

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
    executors=[
        HighThroughputExecutor(
            address=address_by_hostname(),
            label="htex",
            max_workers=1,
            prefetch_capacity=2,
            provider=CobaltProvider(
                queue='CVD_Research',
                account='CVD_Research',
                launcher=AprunLauncher(overrides="-d 64 --cc depth"),
                walltime='3:00:00',
                nodes_per_block=64,
                init_blocks=1,
                min_blocks=1,
                max_blocks=4,
                scheduler_options='#COBALT --attrs enable_ssh=1',
                worker_init='''
module load miniconda-3
source activate /home/lward/exalearn/covid/toxicity-prediction/admet/env
export KMP_AFFINITY=disabled
which python
                ''',
                cmd_timeout=120,
            ),
        ),
    ],
    strategy=None,
)

