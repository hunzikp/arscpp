context("MPLN Sampler")
library(mvtnorm)

set.seed(0)
N <- 100
mu <- rep(3, N)
L <- matrix(rnorm(N*N, 0, 0.01), N, N)
L <- 0.5*(L+t(L)) # L now symmetric
Sigma <- L + diag(N) # Sigma now PD
H <- solve(Sigma)

x <- as.vector(rmvnorm(1, mu, Sigma))
y <- rpois(N, exp(x))

test_that("sampling distribution is correct", {
  set.seed(0)
  mpln <- mplnposterior(y, mu, H)
  smp <- mpln$sample(100, x)
  expect_equal(dim(smp), c(100, N))
  E_x <- colMeans(smp)
  expect_equal(cor(E_x, x) > 0.9, TRUE)
})

test_that("sampling distribution is correct if there are missings", {
  set.seed(0)
  y[sample(length(y), size = 5, replace = FALSE)] <- NA
  mpln <- mplnposterior(y, mu, H)
  smp <- mpln$sample(100, x)
  expect_equal(dim(smp), c(100, N))
  E_x <- colMeans(smp)
  expect_equal(cor(E_x, x) > 0.9, TRUE)
})

test_that("fails if arguments are of unequal size", {
  expect_error(mpln <- plnposterior(1, mu, H))
  expect_error(mpln <- plnposterior(y, 1, H))
  expect_error(mpln <- plnposterior(y, mu, matrix(0, 2, 2)))
})

test_that("fails if starting parameters are of incorrect size", {
  mpln <- mplnposterior(y, mu, H)
  expect_error(smp <- mpln$sample(100, 1))
})

