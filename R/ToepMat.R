MatK = function(x, y ,l, nu){
   ifelse(abs(x - y) > 0, (sqrt(2 * nu) * abs(x - y) / l)^nu /
             (2^(nu - 1) * gamma(nu)) * besselK(x = abs(x - y) * sqrt(2 * nu)/l, nu = nu),
          1.0)
}
