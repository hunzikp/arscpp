#ifndef __LOGDENS_INCLUDED__
#define __LOGDENS_INCLUDED__

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class LogDensity {
  // Abstract base class for Log Densities
public:
  virtual NumericVector h(NumericVector x) = 0;
  virtual NumericVector h_prime(NumericVector x) = 0;
  virtual double h(double x) = 0;
  virtual double h_prime(double x) = 0;
};

class RLogDensity : public LogDensity {
  // Log Density class if functions are passed from R
public:
  Function R_h;
  Function R_h_prime;

  RLogDensity(Function R_h, Function R_h_prime);
  NumericVector h(NumericVector x);
  NumericVector h_prime(NumericVector x);
  double h(double x);
  double h_prime(double x);
};

class PoissonLogNormal : public LogDensity {
  // Log Density of poisson-log-normal posterior
  // log(pois(y; exp(z))) + log(norm(z; mu, sigma^2))
public:
  double y;
  double mu;
  double sigma;

  PoissonLogNormal(double y, double mu, double sigma);

  NumericVector clean(NumericVector x);

  NumericVector h(NumericVector x);
  NumericVector h_prime(NumericVector x);
  double h(double x);
  double h_prime(double x);
};

class NegBinLogNormal : public LogDensity {
  // Log Density of negative-binomial-log-normal posterior
  // log(negbin(y; exp(z), r)) + log(norm(z; mu, sigma^2))
public:
  double y;
  double mu;
  double sigma;
  double r;

  NegBinLogNormal(double y, double mu, double sigma, double r);

  NumericVector clean(NumericVector x);

  NumericVector h(NumericVector x);
  NumericVector h_prime(NumericVector x);
  double h(double x);
  double h_prime(double x);
};

#endif
