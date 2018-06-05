#' Construct an Adaptive Rejection Sampler
#'
#' \code{ars} returns an object for sampling from an arbitrary univariate log-concave density.
#'
#' \code{ars} instantiates a C++ object for sampling from an arbitrary univariate log-concave density using
#' the Adaptive Rejection Sampler by Wild and Gilks (1993, Algorithm AS 287).
#' See the examples for a typical use case.
#'
#' @param f An R function computing the (potentially unnormalized) log-density of the distribution we want to sample from.
#' Must take a numeric vector as the first argument and must be vectorized. Additional arguments to \code{f} may be passed
#' through \code{...}.
#' @param f_prime The first derivative of \code{f} with regard to the variable we want to sample.
#' Must take a numeric vector as the first argument and must be vectorized. Additional arguments to \code{f_prime} may be passed
#' through \code{...}.
#' @param x At least two starting points; at least one starting point should be on the left/right of the mode of the target
#' distribution, respectively. If this is not the case, \code{ars} will try to find valid starting points.
#' @param xlb A lower bound for the support of the target distribution. May be \code{-Inf}.
#' @param xrb An upper bound for the support of the target distribution. May be \code{Inf}.
#' @param max_points Maximum number of support points used to construct the ARS upper/lower hulls. Modifications from the
#' default value should not be necessary.
#' @param ... Additonal named arguments passed to \code{f} and \code{f_prime}.
#'
#' @return A \code{ARS_Rwrapper} object.
#'
#' @examples
#' # Normal log pdf
#' f <- function(x, mu, sd) {
#'   -0.5*(sd^(-2))*(x-mu)^2
#' }
#'
#' # Derivative of normal log pdf
#' f_prime <- function(x, mu, sd) {
#'   -1*(sd^(-2))*(x-mu)
#' }
#'
#' # Make sampler & sample 100 points from standard normal
#' ars_sampler <- ars(f, f_prime, c(-1,1), mu = 0, sd = 1)
#' smp <- ars_sampler$sample(100)
#'
ars <- function(f, f_prime, x, xlb = -Inf, xrb = Inf, max_points = 100, ...) {
  h <- function(x) {f(x, ...)}
  h_prime <- function(x) {f_prime(x, ...)}
  ARS_Rwrapper$new(h, h_prime, x, xlb, xrb, max_points)
}

#' Sample from a Poisson-Log-Normal posterior
#'
#' \code{plnposterior} returns an object for sampling from the PLN posterior.
#'
#' The returned object samples from
#'
#' \deqn{x \sim p(x | y, \mu, sd) \propto dPois(y; exp(x)) \times dNorm(x; \mu, sd)}
#'
#' @param y A numeric vector of counts.
#' @param mu A numeric vector of prior means (must be of same length as \code{y}).
#' @param sd A numeric vector of prior sds (must be of same length as \code{y}).
#'
#' @return A \code{PLNPosterior} object.
#'
#' @examples
#' N <- 100
#' mu <- rep(0, N)
#' sd <- rep(1, N)
#' x <- rnorm(N, mu, sd)
#' y <- rpois(N, exp(x))
#' pln <- plnposterior(y, mu, sd)
#' smp <- pln$sample(100)
#'
plnposterior <- function(y, mu, sd) {
  PLNPosterior$new(y, mu, sd)
}

#' Sample from a multivariate Poisson-Log-Normal posterior
#'
#' \code{mplnposterior} returns an object for sampling from the MPLN posterior
#' using an ARS-in-Gibbs sampler.
#'
#' The returned object samples from
#'
#' \deqn{x \sim p(x | y, \mu, \Sigma) \propto \prod_i [dPois(y_i; exp(x_i))] \times dMVNorm(x; \mu, \Sigma)}
#'
#' @param y A numeric vector of counts of size N.
#' @param mu The prior mean vector of size N.
#' @param H A sparse precision matrix of size NxN, preferably in 'dsCMatrix' (sparse column) format.
#'
#' @return A \code{MPLNPosterior} object.
#'
#' @examples
#' N <- 100
#' mu <- rep(0, N)
#' sd <- rep(2, N)
#' x <- rnorm(N, mu, sd)
#' y <- rpois(N, exp(x))
#' H <- as(diag(1/(sd^2)), 'dsCMatrix')
#' mpln <- mplnposterior(y, mu, H)
#' smp <- mpln$sample(100, x) # x are starting values
#'
mplnposterior <- function(y, mu, H) {
  H <- as(H, 'dgCMatrix')
  MPLNPosterior$new(y, mu, H)
}


