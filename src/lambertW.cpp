/* lambertW.cpp

Original code by Avraham Adler in lamW R package (see below). This version edited by Philipp Hunziker.

Copyright (C) 2015, Avraham Adler
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

References:

Corless, R. M.; Gonnet, G. H.; Hare, D. E.; Jeffrey, D. J. & Knuth, D. E. "On the Lambert W function", Advances in Computational Mathematics,
Springer, 1996, 5, 329-359

Veberič, Darko. "Lambert W function for applications in physics."
Computer Physics Communications 183(12), 2012, 2622-2628

Veberič used for Fritsch iteration step; access to original paper currently unavailable
Need to retain Halley step for -7e-3 < x 7e-3 where the Fritsch may underflow and return NaN
*/

#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

const double EPS = 2.2204460492503131e-16;
const double M_1_E = 1.0 / M_E;

/* Fritsch Iteration as found in http://arxiv.org/pdf/1209.0735.pdf
* W_{n+1} = W_n * (1 + e_n)
* e_n = z_n / (1 + W_n) * (q_n - z_n) / (q_n - 2 * z_n)
* z_n = ln(x / W_n) - W_n
* q_n = 2 * (1 + W_n) * (1 + W_n + 2 / 3 * z_n)
*/

double FritschIter(double x, double w_guess){
  double w = w_guess;
  int MaxEval = 6;
  bool CONVERGED = false;
  double k = 2.0 / 3.0;
  int i = 0;
  do {
    double z = log(x / w) - w;
    double w1 = w + 1.0;
    double q = 2 * w1 * (w1 + k * z);
    double qz = q - z;
    double e = z / w1 * qz / (qz - z);
    CONVERGED = fabs(e) <= EPS;
    w *= (1 + e);
    ++i;
  } while (!CONVERGED && i < MaxEval);
  return(w);
}

/* Halley Iteration
Given x, we want to find W such that Wexp(W) = x, so Wexp(W) - x = 0.
We can use Halley iteration to find this root; to do so it needs first and second derivative.
f(W)    = W * exp(W) - x
f'(W)   = W * exp(W) + exp(W)       = exp(W) * (W + 1)
f''(W)  = exp(W) + (W + 1) * exp(W) = exp(W) * (W + 2)
Halley Step:
W_{n+1} = W_n - {2 * f(W_n) * f'(W_n)} / {2 * [f'(W_n)]^2 - f(W_n) * f''(W_n)}
*/

double HalleyIter(double x, double w_guess){
  double w = w_guess;
  int MaxEval = 12;
  bool CONVERGED = false;
  int i = 0;
  do {
    double ew = exp(w);
    double w1 = w + 1.0;
    double f0 = w * ew - x;
    f0 /= ((ew * w1) - (((w1 + 1.0) * f0) / (2 * w1))); /* Corliss et al. 5.9 */
CONVERGED = fabs(f0) <= EPS;
w -= f0;
++i;
  } while (!CONVERGED && i < MaxEval);
  return(w);
}

double lambertW0_CS(double x) {
  double result;
  double w;
  if (x == std::numeric_limits<double>::infinity()) {
    result = std::numeric_limits<double>::infinity();
  } else if (x < -M_1_E) {
    result = std::numeric_limits<double>::quiet_NaN();
  } else if (fabs(x + M_1_E) < 4 * EPS) {
    result = -1.0;
  } else if (x <= M_E - 0.5) {
    /* Use expansion in Corliss 4.22 to create (3, 2) Pade approximant
    Numerator: -10189 / 303840 * p ^ 3 + 40529 / 303840 * p ^ 2 + 489 / 844 * p - 1
    Denominator: -14009 / 303840 * p^2 + 355 / 844 * p + 1
    Converted to digits to reduce needed operations
    */
    double p = sqrt(2 * (M_E * x + 1));
    double Numer = ((-0.03353409689310163 * p + 0.1333892838335966) * p + 0.5793838862559242) * p - 1;
    double Denom = (-0.04610650342285413 * p + 0.4206161137440758) * p + 1;
    w = Numer / Denom;
    if (fabs(x) <= 7e-3) {
      /* Use Halley step near 0 as this version of Fritsch may underflow */
      result = HalleyIter(x, w);
    } else {
      result = FritschIter(x, w);
    }
  } else {
    /* Use first five terms of Corliss et al. 4.19 */
    w = log(x);
    double L_2 = log(w);
    double L_3 = L_2 / w;
    double L_3_sq = L_3 * L_3;
    w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq / w + L_3_sq * L_3 / 3;
    result = FritschIter(x, w);
  }
  return(result);
}

double lambertWm1_CS(double x){
  double result;
  double w;
  if (x == 0) {
    result = -std::numeric_limits<double>::infinity();
  } else if (x < -M_1_E || x > 0.0) {
    result = std::numeric_limits<double>::quiet_NaN();
  } else if (fabs(x + M_1_E) < 4 * EPS) {
    result = -1.0;
  } else {
    /* Use first five terms of Corliss et al. 4.19 */
    w = log(-x);
    double L_2 = log(-w);
    double L_3 = L_2 / w;
    double L_3_sq = L_3 * L_3;
    w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq / w + L_3_sq * L_3 / 3;
    result = FritschIter(x, w);
  }
  return(result);
}
