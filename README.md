# Code-for-Perturbation-Robust-Predictive-Modeling-of-Social-Effects-by-Network-Subspace-GLM
This folder contains code for the paper "Subspace Generalized Linear Model for Noisy Network-linked Data".

## Case Study
The dataset comes from [Changing climates of conflict: A social network experiment in 56 schools](https://www.pnas.org/doi/10.1073/pnas.1514483113). 

### case\_study\_data\_ORG\_conflict
The processed data for analysis, contains information about 8685 students. 

### model\_fitting
Generate model fitting results in Table 7-9 and Figure 1, Figure 2, Figure 3, Figure 5, Figure 6, Figure 7.

### N2V\_model\_fitting
Generate model fitting results with network embedding.

### CV\_plot and CV\_case\_study
Generate model validation results presented in Figure 3, Figure 4.

### CV\_plot\_embedding and CV\_case\_study\_embed
Generate model validation results for Node2vec/Diff2vec/Deepwalk presented in Figure 6.

## Simulation
The simulation folder contains two part corresponds to two subsection in our paper. It is worth mentioning that:

1. The dimension r=1 of intersection subspace is assumed to be known for our simulation

2. We adjust the signals before the eigenvectors for prediction so that they are consistent for networks under the same setting but with different density.

3. The python code to generate network embedding results with diff2vec from package karateclub is required to be running under python version 3.8-3.9.

There are many .R file to generate simulation results under different setting. For example:

### SBM\_SP
Generate the simulation results of our method in Table 1-4 under SBM. 


### subspace\_SBM\_deepwalk
Generate embedding results of deepwalk method under SBM.

### embedding\_Logit
Generate the simulation results in Table 5 with network embedding results above.
