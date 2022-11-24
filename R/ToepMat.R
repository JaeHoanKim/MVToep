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

rmvToep = function(n, mu = NULL, Sigma, eigtol = 1e-8, symtol = 1e-8){
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
   lambda = eigen(Sigma)$values
   if (min(lambda) < eigtol){
      stop("Sigma is ill-conditioned or not Positive definite; try rmvMat or rmvRBF if applicable.")
   } else{
      out = matrix(0, n, N)
      for (k in 1:n){
         vec = rep(0, N)
         vec[1] = sqrt(lambda[1]) * rnorm(1) / sqrt(N)
         vec[(N / 2) + 1] = sqrt(lambda[(N / 2) + 1]) * rnorm(1) / sqrt(N)
         i=sqrt(as.complex(-1))
         for(j in 2:(N/2)){
            uj = rnorm(1); vj = rnorm(1)
            vec[j] = (sqrt(lambda[j]) * (uj + i * vj)) / (sqrt(2 * N))
            vec[N + 2 - j] = (sqrt(lambda[j]) * (uj - i * vj)) / (sqrt(2 * N))
         }
         out[k, ] = Re(fft(vec))
      }
      out = out + matrix(mu, n, N, byrow = T)
   }
   return(out)
}
