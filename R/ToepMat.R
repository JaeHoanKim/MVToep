#' Matern covariance function
#'
#' @param x vector input
#' @param y vector input (same length with x)
#' @param l length scale (non-negative)
#' @param nu smoothness parameter (non-negative)
#'
#' @return
#' @export
#'
#' @examples
MatK = function(x, y ,l, nu){
   if (length(x) != length(y)){
      stop("the length of x and y should be equal")
   }
   ifelse(abs(x - y) > 0, (sqrt(2 * nu) * abs(x - y) / l)^nu /
             (2^(nu - 1) * gamma(nu)) * besselK(x = abs(x - y) * sqrt(2 * nu)/l, nu = nu),
          1.0)
}
