#' Matern covariance function
#'
#' @param x vector input
#' @param y vector input (same length with x)
#' @param l length scale (non-negative)
#' @param nu smoothness parameter (non-negative)
#'
#' @return calculated kernel function value
#' @export
#'
#' @examples
#' MatK(c(1, 2, 3), c(3, 4, 6), 1, 0.5)
#' MatK(1, 1, 1, 1)
MatK = function(x, y ,l, nu){
   if (length(x) != length(y)){
      stop("the length of x and y should be equal")
   }
   ifelse(abs(x - y) > 0, (sqrt(2 * nu) * abs(x - y) / l)^nu /
             (2^(nu - 1) * gamma(nu)) * besselK(x = abs(x - y) * sqrt(2 * nu)/l, nu = nu),
          1.0)
}

#' Multivariate normal samples from Toeplitz-structured Covariance matrix
#'
#' @param n number of samples extracted from the distribution
#' @param Sigma covarinace matrix of the distribution. It should be a Toeplitz matrix
#' @param mu the mean vector of the distribtuion
#' @param eigtol numerical threshold for Positive definiteness in Sigma
#' @param symtol numerical threshold for checking symmetry in Sigma
#'
#' @return n by N (the number of rows in Sigma) matrix, n multivariate normal vectors stacked vertically
#' @export
#'
#' @examples
#' rmvToep(5, Sigma = diag(10), mu = rep(1, 10))
#' rmvToep(10, Sigma = diag(abs(rnorm(6))), mu = rep(1, 6))
rmvToep = function(n, Sigma, mu = NULL, symtol = 1e-8){
   N = dim(Sigma)[1]
   if (dim(Sigma)[2] != N){
      stop("Sigma should be square matrix!")
   }
   if (max(abs(Sigma - t(Sigma))) > symtol){
      stop("Sigma should be symmetric matrix!")
   }
   if (is.null(mu)){
      mu = rep(0, N)
   }
   # lambda = eigen(Sigma)$values
   Sigma_vec = Sigma[1, ] # first row of Sigma, which contains all the values in Sigma
   C_vec = c(Sigma_vec[1:N], Sigma_vec[N:1])
   lambda = Re(fft(C_vec))
   if (min(lambda) < 0){
      stop("Sigma is ill-conditioned or not Positive definite; try rmvMat or rmvRBF if applicable.")
   } else{
      m = 2*N
      out = matrix(0, n, m)
      for (k in 1:n){
         vec = rep(0, m)
         vec[1] = sqrt(lambda[1]) * rnorm(1) / sqrt(m)
         vec[(m / 2) + 1] = sqrt(lambda[(m / 2) + 1]) * rnorm(1) / sqrt(m)
         i=sqrt(as.complex(-1))
         for(j in 2:(m/2)){
            uj = rnorm(1); vj = rnorm(1)
            vec[j] = (sqrt(lambda[j]) * (uj + i * vj)) / (sqrt(2 * m))
            vec[m + 2 - j] = (sqrt(lambda[j]) * (uj - i * vj)) / (sqrt(2 * m))
         }
         out[k, ] = Re(fft(vec))
      }
      out = out[, 1:N]
      out = out + matrix(mu, n, N, byrow = T)
   }
   return(out)
}
