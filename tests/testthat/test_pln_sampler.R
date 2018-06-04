context("PLN Sampler")

N <- 1000
mu <- rep(3, N)
sd <- rep(1, N)

set.seed(0)
x <- rnorm(N, mu, sd)
y <- rpois(N, exp(x))

test_that("sampling distribution is correct", {
  set.seed(0)
  pln <- plnposterior(y, mu, sd)
  smp <- pln$sample(100)
  expect_equal(dim(smp), c(100, N))
  E_x <- colMeans(smp)
  expect_equal(cor(E_x, x) > 0.9, TRUE)
})

test_that("fails if arguments are of unequal size", {
  expect_error(pln <- plnposterior(1, mu, sd))
  expect_error(pln <- plnposterior(y, 1, sd))
  expect_error(pln <- plnposterior(y, mu, 1))
})
