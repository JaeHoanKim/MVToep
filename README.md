
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast Sampling from Multivariate Normal distribution with Toeplitz-structured covariance matrix (MVToep)

<!-- badges: start -->

[![R-CMD-check](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JaeHoanKim/MVToep/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/JaeHoanKim/MVToep/branch/master/graph/badge.svg)](https://codecov.io/gh/JaeHoanKim/MVToep?branch=master)

<!-- badges: end -->

## Intended use

Multivariate normal (Gaussian) sampling from the Toeplitz-structured
covariance matrix is in diverse fields: from obtaining the approximate
solution of the stochastic partial differential equation (see \[1\]) to
directly testing the performance in the simulation studies (see \[2\]).
Therefore, this package intends to efficiently draw samples from the
Toeplitz-structured covariance matrix, with several user-friendly
functions. For handy usage, it provides simple implementations for
sampling from the regular grid points with radial basis function (RBF)
and Matern kernels, with the covariance matrix of the autoregressive
model. The main idea for this algorithm is suggested on \[3\].

## Installation

You can install the development version of MVToep package by:

``` r
# install.packages("devtools")
devtools::install_github("JaeHoanKim/MVToep")
```

Then, use

    library(MVToep)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MVToep)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:

## References

-   Lindgren, Finn, Håvard Rue, and Johan Lindström. “An explicit link
    between Gaussian fields and Gaussian Markov random fields: the
    stochastic partial differential equation approach.” Journal of the
    Royal Statistical Society: Series B (Statistical Methodology) 73.4
    (2011): 423-498.

-   Kim, Jaehoan, Hoyoung Park, and Junyong Park. “High dimensional
    discriminant rules with shrinkage estimators of covariance matrix
    and mean vector.” arXiv preprint arXiv:2211.15063 (2022).

-   Wood, Andrew TA, and Grace Chan. “Simulation of stationary Gaussian
    processes in \[0, 1\]^d.” Journal of computational and graphical
    statistics 3.4 (1994): 409-432.
