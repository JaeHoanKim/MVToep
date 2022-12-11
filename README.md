<!-- badges: start -->

[![R-CMD-check](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/JaeHoanKim/MVToep/branch/master/graph/badge.svg)](https://codecov.io/gh/JaeHoanKim/MVToep?branch=master)

<!-- badges: end -->

## Fast Sampling from Multivariate Normal distribution with Toeplitz-structured covariance matrix (MVToep)

### Intended use

Multivariate normal (Gaussian) sampling from the Toeplitz-structured covariance matrix is in diverse fields: from obtaining the approximate solution of the stochastic partial differential equation (see [1]) to directly testing the performance in the simulation studies (see [2]). Therefore, this package intends to efficiently draw samples from the Toeplitz-structured covariance matrix, with several user-friendly functions. For handy usage, it provides simple implementations for sampling from the regular grid points with radial basis function (RBF) and Matern kernels, with the covariance matrix of the autoregressive model. The main idea for this algorithm is suggested on [3].

### Installation instructions

For the installation of the MVToep package, use the github install function. Type

``` install1
devtools::install_github("JaeHoanKim/MVToep")
``` 

into the R console. Then, use 

``` install2
library(MVToep)
``` 

to use the functions in MVToep.

### Examples

Now that rmvToep function is available, rmvRBF and rmvMat, which are user-friendly functions for regular grids with corresponding covariance kernels would be constructed. Out of developing codes, formal citation of the paper (Wood and Chan; 1994) and making vignette is needed.


### References

- Lindgren, Finn, Håvard Rue, and Johan Lindström. "An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 73.4 (2011): 423-498.

- Kim, Jaehoan, Hoyoung Park, and Junyong Park. "High dimensional discriminant rules with shrinkage estimators of covariance matrix and mean vector." arXiv preprint arXiv:2211.15063 (2022).

- Wood, Andrew TA, and Grace Chan. "Simulation of stationary Gaussian processes in [0, 1]^d." Journal of computational and graphical statistics 3.4 (1994): 409-432.
