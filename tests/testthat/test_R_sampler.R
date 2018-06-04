context("Generic ARS Sampler")

# Normal log pdf
f <- function(x, mu, sd) {
  -0.5*(sd^(-2))*(x-mu)^2
}

# Derivative of normal log pdf
f_prime <- function(x, mu, sd) {
  -1*(sd^(-2))*(x-mu)
}

test_that("generic sampling distribution is correct", {
  ars_sampler <- ars(f, f_prime, c(-1,1), mu = 0, sd = 1)
  set.seed(0)
  smp <- ars_sampler$sample(100)
  expect_equal(shapiro.test(smp)$p.value > 0.05, TRUE)  # Test normality
})

test_that("sampling distribution with lower bound works", {
  ars_sampler <- ars(f, f_prime, c(-1,1), xlb = -2, mu = 0, sd = 1)
  set.seed(0)
  smp <- ars_sampler$sample(100)
  expect_equal(all(smp>-2), TRUE)
})

test_that("sampling distribution with upper bound works", {
  ars_sampler <- ars(f, f_prime, c(-1,1), xrb = 2, mu = 0, sd = 1)
  set.seed(0)
  smp <- ars_sampler$sample(100)
  expect_equal(all(smp<2), TRUE)
})

test_that("automatic determination of lower bound works", {
  ars_sampler <- ars(f, f_prime, c(0.5,1), mu = 0, sd = 1) # Incorrect lower bound
  set.seed(0)
  smp <- ars_sampler$sample(100)
  expect_equal(shapiro.test(smp)$p.value > 0.05, TRUE)  # Test normality
})

test_that("automatic determination of upper bound works", {
  ars_sampler <- ars(f, f_prime, c(-1, -0.5), mu = 0, sd = 1) # Incorrect upper bound
  set.seed(0)
  smp <- ars_sampler$sample(100)
  expect_equal(shapiro.test(smp)$p.value > 0.05, TRUE)  # Test normality
})

test_that("sampler works with extreme starting points", {
  ars_sampler <- ars(f, f_prime, c(-1e7, 1e7), mu = 0, sd = 1)
  set.seed(0)
  smp <- ars_sampler$sample(100)
  expect_equal(shapiro.test(smp)$p.value > 0.05, TRUE)  # Test normality
})

test_that("sampler works with many unordered starting points", {
  set.seed(0)
  ars_sampler <- ars(f, f_prime, rnorm(100), mu = 0, sd = 1)
  set.seed(0)
  smp <- ars_sampler$sample(100)
  expect_equal(shapiro.test(smp)$p.value > 0.05, TRUE)  # Test normality
})

test_that("sampler fails with only one starting point", {
  expect_error(ars_sampler <- ars(f, f_prime, 1, mu = 0, sd = 1))
})

test_that("sampler fails with invalid bounds", {
  expect_error(ars_sampler <- ars(f, f_prime, c(-1, 1), xlb = 0, xrb = -Inf))
})


