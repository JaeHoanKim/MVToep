#' Matern kernel function
#'
#' @param x scalar of vector input
#' @param y scalar vector input (same length with x)
#' @param rho length scale (non-negative)
#' @param nu smoothness parameter (non-negative)
#'
#' @return It returns the calculated Matern kernel function value, along the following formula : \cr
#' \eqn{C(d) = \sigma^2\frac{2^{1 - \nu}}{\Gamma(\nu)} (\sqrt{2\nu}\frac{d}{\rho})^{\nu} K_{\nu}(\sqrt{2\nu}\frac{d}{\rho})}, \cr
#' where \eqn{d = |x-y|}.
#' @export
#'
#' @examples
#' MatK(c(1, 2, 3), c(3, 4, 6), 1, 0.5)
#' MatK(1, 1, 1, 1)
MatK = function(x, y ,rho, nu){
   ifelse(abs(x - y) > 0, (sqrt(2 * nu) * abs(x - y) / rho)^nu /
             (2^(nu - 1) * gamma(nu)) * besselK(x = abs(x - y) * sqrt(2 * nu)/rho, nu = nu),
          1.0)
}

#' RBF kernel function
#'
#' @param l length parameter in the RBF kernel
#' @inheritParams MatK
#'
#' @return It returns the calculated RBF kernel function value, along the following formula : \cr
#' \eqn{C(d) = \exp(-\frac{d^2}{2l^2})}, \cr
#' where \eqn{d = |x-y|}.
#' @export
#'
#' @examples
#' RBFK(c(0, 1, 2), 0, 5)
RBFK = function(x, y, l){
   return(exp(-(x - y)^2 / (2 * l^2)))
}

Sigma.AR.order = function(rho, p, order = NULL){
   if (order = NULL){
      order = p
   }
   R = diag(p)
   for(i in 1:p){
      for(j in 1:p){
         R[i, j] = ifelse(abs(i-j) > order, 0, rho^(abs(i-j)))
      }
   }
   return(R)
}

nnd.C.Toep = function(Sigma){
   N = nrow(Sigma)
   Sigma_vec = Sigma[1, ]
   out = min(Re(fft(c(Sigma_vec[1:N], Sigma_vec[N:1]))))
   return(ifelse(out >= 0, "rmvToep is applicable for the given Sigma!",
                 "min(lambda) < 0; rmvToep cannot be applied for the given Sigma!"))
}

#' Multivariate normal samples from Toeplitz-structured Covariance matrix
#'
#' @param n number of samples extracted from the distribution
#' @param Sigma covarinace matrix of the distribution. It should be a Toeplitz matrix
#' @param mu the mean vector of the distribtuion
#' @param symtol numerical threshold for checking symmetry in Sigma
#'
#' @return n by N (the number of rows in Sigma) matrix, n multivariate normal vectors stacked vertically
#' @export
#'
#' @examples
#' rmvToep(5, Sigma = diag(10), mu = rep(1, 10))
#' rmvToep(10, Sigma = diag(abs(rnorm(6))), mu = rep(1, 6))
rmvToep = function(n, Sigma, mu = rep(0, nrow(Sigma)), symtol = 1e-8){
   N = dim(Sigma)[1]
   if (dim(Sigma)[2] != N){
      stop("Sigma should be a square matrix!")
   }
   if (max(abs(Sigma - t(Sigma))) > symtol){
      stop("Sigma should be a symmetric matrix!")
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

grid_regular_check = function(gridpoints){
   N = length(gridpoints)
   gridpoints = sort(gridpoints)
   point_ini = gridpoints[1]
   gap = gridpoints[2] - gridpoints[1]
   if (max(abs(gridpoints - point_ini - gap * c(0:(N-1)))) > 1e-8){
      stop("gridpoints should be regular!")
   }
}

min_m = function(gridpoints){
   N = length(gridpoints)
   g = ceiling(log(2 * N, 2))
   m = 2^g
   return("m" = m)
}

embed_vec = function(gridpoints, m){
   N = length(gridpoints)
   gridpoints = sort(gridpoints)
   point_ini = gridpoints[1]
   gap = gridpoints[2] - gridpoints[1]
   cj = vector(length = m)
   for (j in 1:m){
      if (j <= (m/2)){
         cj[j] = (j - 1) * gap
      }else{
         cj[j] = (m - (j - 1)) * gap
      }
   }
   return(cj)
}

C.eigval.Mat = function(gridpoints, m, rho, nu){
   vec = embed_vec(gridpoints, m)
   cj = MatK(vec, 0, rho, nu)
   eigval = Re(fft(cj))
   return(list("cj" = cj, "eigval" = eigval))
}

C.eigval.RBF = function(gridpoints, m, l){
   vec = embed_vec(gridpoints, m)
   cj = RBFK(vec, 0, l)
   eigval = Re(fft(cj))
   return(list("cj" = cj, "eigval" = eigval))
}

nnd.C.Mat = function(gridpoints, m, rho, nu){
   out = C.eigval.Mat(gridpoints, m, rho, nu)
   cj = out$cj
   eigval = out$eigval
   if (min(eigval) > 0){
      return(list("cj" = cj, "m" = m, "eigval" = eigval))
   }
   else{
      m = 2 * m
      if (m > 2000){
         stop("It seems that adequate circular matrix is not found. Try smaller l.")
      }
      nnd.C.Mat(gridpoints, m, rho, nu)
   }
}

nnd.C.RBF = function(gridpoints, m, l){
   out = C.eigval.RBF(gridpoints, m, l)
   cj = out$cj
   eigval = out$eigval
   if (min(eigval) > 0){
      return(list("cj" = cj, "m" = m, "eigval" = eigval))
   }
   else{
      m = 2 * m
      if (m > 2000){
         stop("It seems that adequate circular matrix is not found. Try smaller l.")
      }
      nnd.C.RBF(gridpoints, m, l)
   }
}

#' Multivariate normal sampling from RBF covarinace kernel with regular grids
#'
#' @param rho length parameter in Matern kernel
#' @param nu smoothness parameter in Matern kernel
#' @inheritParams rmvRBF
#'
#' @return asdf
#' @export
#'
#' @examples
#' rmvMat(10, c(0:20)/20, 1, 0.5, tau = 1)
rmvMat = function(n, gridpoints, rho, nu, mu = rep(0, length(gridpoints)), tau = 1){
   N = length(gridpoints)
   grid_regular_check(gridpoints)
   embed_result = nnd.C.Mat(gridpoints, m = min_m(gridpoints), rho, nu)
   cj = embed_result$cj
   m = embed_result$m
   lambda = embed_result$eigval
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
   out = out[, 1:N] * tau
   out = out + matrix(mu, n, N, byrow = T)
   return(out)
}

#' Multivariate normal sampling from RBF covarinace kernel with regular grids
#'
#' @param n a
#' @param gridpoints b
#' @param l c
#' @param mu d
#' @param tau e
#'
#' @return f
#' @export
#'
#' @examples
#' rmvRBF(100, c(0:20)/20, l = 0.1, tau = 1)
rmvRBF = function(n, gridpoints, l, mu = rep(0, length(gridpoints)), tau = 1){
   N = length(gridpoints)
   grid_regular_check(gridpoints)
   embed_result = nnd.C.RBF(gridpoints, m = min_m(gridpoints), l)
   cj = embed_result$cj
   m = embed_result$m
   lambda = embed_result$eigval
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
   out = out[, 1:N] * tau
   out = out + matrix(mu, n, N, byrow = T)
   return(out)
}
