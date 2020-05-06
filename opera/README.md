# OPERA 

[OPERA](https://github.com/NIEHS/OPERA) is an open-source tool from the NIH for running chemoinformatics models.
The latest versions of OPERA contain dozens of models including both common physiochemical properties (e.g., melting point)
and, of relevance here, models from the NIH's collaborative efforts on toxicity modeling.
These models include:

- [CERAPP](https://www.ncbi.nlm.nih.gov/pubmed/26908244)/[CoMPRA](https://www.ncbi.nlm.nih.gov/pubmed/32074470): Screening for whether molecules are likely to be endocrine distruptors through whether a molecules interacts with estrogen or androgen receptors.
- [CATMoS](https://ntp.niehs.nih.gov/ntp/about_ntp/sacatm/2019/september/presentations/2-2-mansouri-508.pdf): A combination of several different toxicity measures into a single model

Each model is based on an ensemble of chemoinformatics models produced by a few dozen groups around the world.

## Installation

First, install the conda environment: `conda env create --file environment.yml`.

Then, you will need to install OPERA. Our code is written to use the parallel, CLI version of OPERA available [here](https://github.com/NIEHS/OPERA/releases/tag/v2.5-beta2).

The installer for Opera requires you have Java installed on your system and will install the MATLAB runtime environment if you do not have it.
Installation uses a GUI and will install the OPERA code and a MATLAB runtime in separate directories.
Copy the directories produced from the installation to whatever system you need.

Then, edit the `set_envs.sh` file to point to the install paths of thse folders and 
the location of the `libjvm.so` library on your system.

Finally, download [CDK](https://sourceforge.net/projects/cdk/files/cdk/1.4.15/cdk-1.4.15.jar/download) and update the 
`set_envs.sh` file on your system to point the location of the JAR file.

Your first time running the code will create a `~/.mcrCache*/` directory and then produce an error.
Edit the `OPERA_installdir.txt` file in this directory to contain the path to your OPERA installation.

### Running OPERA on LCFs

There are a few challenges to running OPERA on an LCF. 

First, make sure your environment has an installation of Java and 
make sure your Parsl executor loads the module for each worker.

Then you must edit the `run_OPERA_P.sh` file in your OPERA installation
to ensure that each compute node uses a different `~/.mcrCache*/` directory.
Do so by first running the `parsl_inference.py` on a small dataset on a login node
to fill the `~/.mcrCache/` directory with all of the needed cache files.
Next, add the following lines to your `run_OPERA_P.sh` script:

```bash
# Set the cache location
echo "Changing the cache root"
export MCR_CACHE_ROOT=`pwd`/`hostname`/
mkdir -p $MCR_CACHE_ROOT
cp -r ~/.mcrCache9.4/ $MCR_CACHE_ROOT
```

A complete example of the `run_OPERA_P.sh` showing where to place these lines is available
in [this repo](./docs/run_OPERA_P.sh).

If you need to alter the Parsl configuration for a new system (mine is tuned for Theta),
there are a few stumbling blocks I encountered:

1. The application does not scale to large numbers of cores per node. Use more than one worker per node
   and manually set the number of cores per node to prevent oversubscription.
1. The application uses a lot of memory. I find that I can fit only 4 workers per node with batch
   sizes of 1024 molecules per task.

