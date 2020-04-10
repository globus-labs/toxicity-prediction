# Tox21 Screening Results

This folder contains CSV files that are the result of screening based on a model trained on the Tox21 dataset.
The Tox21 dataset was produced in 2014 by the NIH and contains the output of 12 different assays, 
8 of which deal with "Nuclear Receptor Signaling" and 4 of which deal with stress response. 
The full details and links to descriptions of the assays can be found [on the NIH's website](https://tripod.nih.gov/tox21/challenge/data.jsp).

We used a graph-convolution model trained via DeepChem, 
which was the best-performing model on a [benchmarking effort
performed by the DeepChem developers](http://dx.doi.org/10.1039/C7SC02664A).

The first 12 columns of the CSV file are the screening outputs
where each column is the probability of a molecular being 
found as "active" by a certain assay. 
Interpret values close to 1 a liklihood of the molecule being "toxic."
The 13th column is the SMILES string of the molecule which was screened.

The remaining columns are debugging and performance information.
They include at least the FuncX task ID associated with a batch of results, 
the time required to evaluate the batch,
start and end times of execution,
the Theta node which executed it,
and the number of cores used on that node.
