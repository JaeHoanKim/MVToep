test_that("MatK", {
   expect_equal(MatK(1, 1, 1, 1), 1)
   expect_equal(MatK(c(1, 2, 3), 1, rho = 2, nu = 3), c(1.0000000, 0.83910663, 0.53592547))
})

test_that("RBFK", {
   expect_equal(RBFK(1, 1, 1), 1)
   expect_equal(RBFK(c(1, 2, 3), 1, l = 3), c(1.0000000, 0.945959469, 0.800737403))
})

## check the size of the output

Sigma1 = matrix(c(1, 0.1, 0, 0, 0,
                  0.1, 1, 0.1, 0, 0,
                  0, 0.1, 1, 0.1, 0,
                  0, 0, 0.1, 1, 0.1,
                  0, 0, 0, 0.1, 1), 5, 5)

test_that("rmv dimension test", {
   expect_equal(dim(rmvToep(50, Sigma1, mu = c(1, 2, 3, 4, 5))), c(50, 5))
   expect_equal(dim(rmvMat(30, c(0:50)/50, rho = 1, nu = 1)[[1]]), c(30, 51))
   expect_equal(dim(rmvMat(1, c(0:50)/50, rho = 1, nu = 1)[[1]]), c(1, 51))
   expect_equal(dim(rmvRBF(50, c(0:10), l = 0.1, mu = 10 * runif(10))[[1]]), c(50, 11))
})

## check the result for AR matrix and nonnegativity check

Sigma11 = Sigma.AR.order(0.3, 20)
Sigma12 = Sigma.AR.order(0.8, 10, order = 2)

test_that("rmvToep result test", {
   expect_equal(dim(Sigma11), c(20, 20))
   expect_equal(nnd.C.Toep(Sigma11), "rmvToep is applicable for the given Sigma!")
   expect_equal(nnd.C.Toep(Sigma12), "min(lambda) < 0; rmvToep cannot be applied for the given Sigma!")
})

## error expected for ill conditioned matrix

test_that("error check for ill-conditioned matrix",{
   expect_error(rmvRBF(10, c(0:30)/30, l=3), "It seems that adequate circular matrix is not found. Try smaller l.")
   expect_error(rmvMat(10, c(0:50)/50, rho = 3, nu = 1), "It seems that adequate circular matrix is not found. Try smaller l.")
})

## exact error location check

test_that("error in nnd function", {
   expect_error(nnd.C.RBF(c(0:30)/30, m = 64, l = 3), "It seems that adequate circular matrix is not found. Try smaller l.")
   expect_error(nnd.C.Mat(c(0:50)/50, m = 128, rho = 3, nu = 1), "It seems that adequate circular matrix is not found. Try smaller l.")
})

## compatibility check

Sigma2 = matrix(c(2, 1, 0,
                 1, 2, 1,
                 0, 1.001, 2,
                 1, 2, 3), 3, 4)

Sigma3 = matrix(c(2, 1, 0,
                  1, 2, 1,
                  0, 1.001, 2), 3, 3)

Sigma4 = matrix(c(1, 0.99, 0.97, 0.94, 0.89,
                  0.99, 1, 0.99, 0.97, 0.94,
                  0.97, 0.99, 1, 0.99, 0.97,
                  0.94, 0.97, 0.99, 1, 0.99,
                  0.89, 0.94, 0.97, 0.99, 1), 5, 5)
vec1 = abs(rnorm(25))
Sigma5 = matrix(vec1, 5, 5) + t(matrix(vec1, 5, 5))


grid1 = c(0:8)/8
grid1[3] = grid1[3] + 0.001

test_that("Compatibility check", {
   expect_error(rmvToep(5, Sigma2), "Sigma should be a square matrix!")
   expect_error(rmvToep(5, Sigma3), "Sigma should be a symmetric matrix!")
   expect_error(rmvToep(5, Sigma5), "Sigma is a symmetric matrix, but not a Toeplitz matrix!")
   expect_error(rmvToep(5, Sigma4), "Sigma is ill-conditioned or not Positive definite; try rmvMat or rmvRBF if applicable.")
   expect_error(grid_regular_check(grid1), "gridpoints should be regular!")
   expect_warning(grid_regular_check(c(1, 2, 3, 5, 4)), "The gridpoints are not sorted. The gridpoints will be used in an ascending order!")
   expect_error(Sigma.AR.order(0.3, 20, -1), "order should be given as a non-negative integer!")
})
