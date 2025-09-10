# spatialPCA
This repository contains R-code related to the manuscript "New M-estimator of the leading principal component" by Virta, Radojičić and Voutilainen (2025+).

The manuscript investigates the minimization of the objective function $v \mapsto \mathrm{E}\| X - v \| \| X + v \|$ and shows that the minimizers $v_0$ can be interpreted as estimating the direction of the leading principal component of $X$.

Preprint available in arXiv (add a link)

All code is in R.

The code files are:

- code_algorithm.R, an implementation of the algorithm for computing the mininizer as described in Section 4 of the manuscript

As proven in Corollary 2 of the manuscript, the algorithm is guaranteed to converge in a finite number of iterations.

## Authors

Virta J., Radojičić U. and Voutilainen M.

## License

GNU GPLv3
