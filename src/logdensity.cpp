#include <Rcpp.h>
#include "logdensity.h"
using namespace Rcpp;
using namespace std;


// RLogDensity implementation
// -> Log Density class if functions are passed from R
RLogDensity::RLogDensity(Function R_h, Function R_h_prime) :
    R_h(R_h), R_h_prime(R_h_prime) {}

NumericVector RLogDensity::h(NumericVector x) {
  NumericVector rv = R_h(x);
  return rv;
}

NumericVector RLogDensity::h_prime(NumericVector x) {
  NumericVector rv = R_h_prime(x);
  return rv;
}

double RLogDensity::h(double x) {
  NumericVector rv = R_h(x);
  return rv.at(0);
}

double RLogDensity::h_prime(double x) {
  NumericVector rv = R_h_prime(x);
  return rv.at(0);
}


// PoissonLogNormal implementation
PoissonLogNormal::PoissonLogNormal(double y, double mu, double sigma) :
  y(y), mu(mu), sigma(sigma) {}

NumericVector PoissonLogNormal::clean(NumericVector x) {
  // ensure x is below 700 (max value for which we can compute exp)
  for (int i = 0; i < x.size(); ++i) {
    if (x.at(i) > 700) {
      x.at(i) = 700;
    }
  }
  return x;
}

NumericVector PoissonLogNormal::h(NumericVector x) {
  x = clean(x);
  NumericVector rv = x*y - exp(x) - 0.5*pow(sigma, -2)*pow(x - mu, 2);
  return rv;
}

NumericVector PoissonLogNormal::h_prime(NumericVector x) {
  x = clean(x);
  NumericVector rv = y - exp(x) - pow(sigma, -2)*(x - mu);
  return rv;
}

double PoissonLogNormal::h(double x) {
  NumericVector rv = h(wrap(x));
  return rv.at(0);
}

double PoissonLogNormal::h_prime(double x) {
  NumericVector rv = h_prime(wrap(x));
  return rv.at(0);
}


