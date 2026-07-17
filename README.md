# spatialPCA
This repository contains R-code related to the manuscript "New M-estimator of the leading principal component" by Virta, Radojičić and Voutilainen (2025+).

The manuscript investigates the minimization of the objective function $v \mapsto \mathrm{E}\| X - v \| \| X + v \|$ and shows that the minimizers $v_0$ can be interpreted as estimating the direction of the leading principal component of $X$.

Preprint available in [arXiv](https://arxiv.org/abs/2510.02799).

All code is in R.

The code files are:

- code_algorithm.R, an implementation of the algorithm for computing the mininizer as described in Section 4 of the manuscript.
- code_simulation_1.R, code for running Simulation study #1 in Section 5.1 of the manuscript and for producing the result Figures 2 and 3.
- code_simulation_2.R, code for running Simulation study #2 in Section 5.2 of the manuscript and for producing the result Figure 4. For reproducing Figure 4 without rerunning the simulation, use "res_cov.csv"
- code_simulation_3.R, code for running the additional experiment in Section 6 and producing the result Figure 5 about estimating PC2 with the method.
- code_simulation_4.R, code for reproducing the sensitivity study within Simulation 2 and producing Table 1. 

As proven in Corollary 2 of the manuscript, the algorithm is guaranteed to converge in a finite number of iterations.

## Authors

Virta J., Radojičić U. and Voutilainen M.

## License

GNU GPLv3
