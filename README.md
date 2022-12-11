<!-- badges: start -->

[![R-CMD-check](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/JaeHoanKim/MVToep/branch/master/graph/badge.svg)](https://codecov.io/gh/JaeHoanKim/MVToep?branch=master)

<!-- badges: end -->

## Fast Sampling from Multivariate Normal distribution with Toeplitz-structured covariance matrix (MVToep)

### Intended use

Multivariate normal (Gaussian) sampling from the Toeplitz-structured covariance matrix is in diverse fields: from obtaining the approximate solution of the stochastic partial differential equation to directly testing the performance in the simulation studies. Therefore, this package intends to efficiently draw samples from the Toeplitz-structured covariance matrix, with several user-friendly functions. For handy usage, it provides simple implementations for sampling from the regular grid points with radial basis function (RBF) and Matern kernels, with the covariance matrix of the autoregressive model. The main idea for this algorithm is suggested on Wood and Chan (1994).

### Installation instructions

For the installation of the MVToep package, use the github install function. Type *devtools::install_github("JaeHoanKim/MVToep")* into the R console. Then, use *library(MVToep)* to use two functions in MVToep.

### Remainder of the project

Now that rmvToep function is available, rmvRBF and rmvMat, which are user-friendly functions for regular grids with corresponding covariance kernels would be constructed. Out of developing codes, formal citation of the paper (Wood and Chan; 1994) and making vignette is needed.
