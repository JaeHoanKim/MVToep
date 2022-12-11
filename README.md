
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

The following is a basic example which shows how to use `rmvToep`
function.

To use the `rmvToep` function, a symmetric and Toeplitz-structured
covariance matrix should be provided. `mvToep` package provides one
class of covariance matrix that satisfies these conditions in
`Sigma.AR.order` function, which is the covariance matrix of the
autoregressive (AR) model. One can control the number of nonzero
elements in this matrix by changing the `order` argument.

``` r
library(MVToep)
Sigma = Sigma.AR.order(0.7, 200, order = 5)
print(Sigma[1:8, 1:8])
#>         [,1]    [,2]    [,3]   [,4]   [,5]    [,6]    [,7]    [,8]
#> [1,] 1.00000 0.70000 0.49000 0.3430 0.2401 0.16807 0.00000 0.00000
#> [2,] 0.70000 1.00000 0.70000 0.4900 0.3430 0.24010 0.16807 0.00000
#> [3,] 0.49000 0.70000 1.00000 0.7000 0.4900 0.34300 0.24010 0.16807
#> [4,] 0.34300 0.49000 0.70000 1.0000 0.7000 0.49000 0.34300 0.24010
#> [5,] 0.24010 0.34300 0.49000 0.7000 1.0000 0.70000 0.49000 0.34300
#> [6,] 0.16807 0.24010 0.34300 0.4900 0.7000 1.00000 0.70000 0.49000
#> [7,] 0.00000 0.16807 0.24010 0.3430 0.4900 0.70000 1.00000 0.70000
#> [8,] 0.00000 0.00000 0.16807 0.2401 0.3430 0.49000 0.70000 1.00000
```

Before applying `rmvToep` function on `Sigma`, one can check the
availability of this algorithm by using `nnd.C.Toep` function. This
function returns the comment whether one can apply `rmvToep` function or
not. `Sigma` should be positive semi-definite matrix to apply `rmvToep`
function.

``` r
nnd.C.Toep(Sigma)
#> [1] "rmvToep is applicable for the given Sigma!"
```

Here is the example in which one cannot apply `rmvToep` function:

``` r
Sigma2 = Sigma.AR.order(0.7, 200, order = 1)
nnd.C.Toep(Sigma2)
#> [1] "min(lambda) < 0; rmvToep cannot be applied for the given Sigma!"
```

The example below illustrates the sampling result from `rmvToep`
function.

``` r
sample1 = rmvToep(500, Sigma)
plot(sample1[, 1], sample1[, 2])
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="50%" style="display: block; margin: auto;" />

``` r
(cor(sample1[, 1], sample1[, 2]))
#> [1] 0.7079458
```

Simply saying, `order` argument indicates the nonzero bandwidth of the
matrix. Under the AR model, the correlation coefficient of the
neighborhood is `rho`, which matches the above result.

``` r
plot(sample1[, 1], sample1[, 21])
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="50%" style="display: block; margin: auto;" />

``` r
(cor(sample1[, 1], sample1[, 21]))
#> [1] -0.03113447
```

In contrast, two indexes located far from each other show negligable
correlation.

`rmvMat` function works in a simpler way. This is the function to draw
multivariate normal samples under Matern covariance kernel at a regular
grid, which can be directly applied to obtain the approximate solution
of the stochastic partial differential equation (SPDE). The object
`result_Mat` contains two arguments as a separate list. The first list
is the multivariate samples, and the second the used covariance matrix.
One can use `?MatK` to check how the covariance matrix is calculated.

``` r
grid = c(0:300)/50
result_Mat = rmvMat(n = 40, grid, rho = 3, nu = 0.5, tau = 2)
sample2 = result_Mat[[1]]
```

By the visualization below, one can indirectly check that each row of
the `sample2` indicates the Gaussian process with Matern covariance
matrix.

``` r
plot(sample2[1, ], ylim = c(-6, 6))
for (i in 1:7){
   lines(sample2[i, ], col = i)
}
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="50%" style="display: block; margin: auto;" />

`rmvRBF` works similarly as `rmvMat` function. The only difference is
the required arguments for the kernel setting. It is worth noting that
rmvRBF is more likely to suffer from the ill-conditionedness of the
covariance matrix due to its nature. One can use `?RBFK` to check how
the covariance matrix is calculated.

``` r
grid = c(0:13)/13
result_RBF = rmvRBF(n = 40, grid, l = 0.1, tau = 1)
sample3 = result_RBF[[1]]
plot(sample3[1, ], ylim = c(-3, 3))
for (i in 1:7){
   lines(sample3[i, ], col = i)
}
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="50%" style="display: block; margin: auto;" />

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
