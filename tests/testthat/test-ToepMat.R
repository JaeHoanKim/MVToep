test_that("MatK", {
   expect_equal(MatK(1, 1, 1, 1), 1)
   expect_equal(MatK(c(1, 2, 3), 1, rho = 2, nu = 3), c(1.0000000, 0.83910663, 0.53592547))
})

test_that("RBFK", {
   expect_equal(RBFK(1, 1, 1), 1)
   expect_equal(RBFK(c(1, 2, 3), 1, l = 3), c(1.0000000, 0.945959469, 0.800737403))
})

## error expected for ill conditioned matrix

test_that("error check for ill-conditioned matrix",{
   expect_error(rmvRBF(10, c(0:30)/30, l=3), "It seems that adequate circular matrix is not found. Try smaller l.")
})
