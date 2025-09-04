## Eigenvalue gap finder

Contains the code for reproducing the experiments in the upcoming arXiv preprint [1]. The main function `gapfinder_main` aims to find gaps in the spectrum of a real symmetric matrix $A$, 
by approximating the traces of spectral projectors associated with sections of its spectrum, using a combination of Hutchinson's stochastic trace estimator and the Lanczos algorithm for computing quadratic forms. 

[1] Michele Benzi, Michele Rinelli, Igor Simunec, *Estimation of spectral gaps for sparse symmetric matrices*, to appear on arXiv, 2024.
