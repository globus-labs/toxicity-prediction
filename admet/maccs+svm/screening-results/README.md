# ADMET Screening Data

This folder contains the screening results from using models
trained on the [admetSAR](http://lmmd.ecust.edu.cn/admetsar2) database.

The models are SVM models that use FP4 and MACCS fingerprints as inputs.

Each CSV file contains the SMILES string for each molecule screened, 
whatever identifier was included with it in the original data file,
some debug information that measures the performance on Theta,
and toxicity predictions.
The toxicity screening columns are probability of each compound being toxic
and correspond to models trained on the following datasets:

- T_hERG_I: [human Ether-à-go-go-Related Gene inhibitor, dataset 1](https://onlinelibrary.wiley.com/doi/10.1002/minf.201000159)
- T_hERG_II: [human Ether-à-go-go-Related Gene inhibitor, dataset 2](https://doi.org/10.1021/mp300023x)
- T_Carc_I: [Chemical carcinogen in rats](http://doi.org/10.1002/qsar.200860192)
- T_AMES_I: [Ames test for mutagenicity](https://doi.org/10.1021/ci300400a)

*TODO*: Get these descriptions reviewed and edited by an expert 
