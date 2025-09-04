## Eigenvalue gap finder

Contains the code for reproducing the experiments in the arXiv preprint [1]. The main function `gapfinder_main` aims to find gaps in the spectrum of a real symmetric matrix $A$, 
by approximating the traces of spectral projectors associated with sections of its spectrum, using a combination of Hutchinson's stochastic trace estimator and the Lanczos algorithm for computing quadratic forms. 
The test matrix used in `test_gapfinder_dft` is generated with `h2chain-matrices.py` using the Python package `pyscf`.

[1] Michele Benzi, Michele Rinelli, Igor Simunec, *Estimation of spectral gaps for sparse symmetric matrices*, arXiv:2410.15349, 2024.
