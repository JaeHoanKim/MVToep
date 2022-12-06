<!-- badges: start -->
  [![R-CMD-check](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Multivariate normal distribution from Toeplitz-structured covariance matrix (MVToep)
------------------------------------------------------------------------

### Intended use

This package is intend to efficiently draw multivariate normal samples from the Toeplitz-structured covariance matrices, which is widely used throughout the fields. For handy usage, it provides simple implementations for sampling from the regular grid points with radial basis function (RBF) and Matern kernels. The main idea for this algorithm is based on Wood and Chan (1994).

### Installation instructions

For the installation of the MVToep package, use the github install function. Type *devtools::install_github("JaeHoanKim/MVToep")* into the R console. Then, use *library(MVToep)* to use two functions in MVToep.

### Remainder of the project

Now that rmvToep function is available, rmvRBF and rmvMat, which are user-friendly functions for regular grids with corresponding covariance kernels would be constructed. Out of developing codes, formal citation of the paper (Wood and Chan; 1994) and making vignette is needed.
