test_that("MatK", {
   expect_equal(MatK(1, 1, 1, 1), 1)
   expect_equal(MatK(c(1, 2, 3), 1, rho = 2, nu = 3), c(1.0000000, 0.83910663, 0.53592547))
})

test_that("RBFK", {
   expect_equal(RBFK(1, 1, 1), 1)
   expect_equal(RBFK(c(1, 2, 3), 1, l = 3), c(1.0000000, 0.945959469, 0.800737403))
})

## check the size of the output

test_that("rmvMat dimension", {
   expect_equal(dim(rmvMat(30, c(0:50)/50, rho = 1, nu = 1)), c(30, 51))
   expect_equal(dim(rmvMat(1, c(0:50)/50, rho = 1, nu = 1)), c(1, 51))
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

Sigma1 = matrix(c(2, 1, 0,
                 1, 2, 1,
                 0, 1.001, 2))

grid1 = c(0:8)/8
grid1[3] = grid1[3] + 0.001

test_that("Compatibility check", {
   expect_error(rmvToep(5, Sigma1), "Sigma should be a square matrix!")
   expect_error(grid_regular_check(grid1), "gridpoints should be regular!")
})



