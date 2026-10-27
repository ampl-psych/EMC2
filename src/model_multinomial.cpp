#include <Rcpp.h>
#include <cmath>
#include "model_multinomial.h"
#include "hcubature.h"

namespace {

struct ProbitChoice {
  const double* utility;
  int n, winner;
};

int probit_integrand(unsigned, const double* x, void* data, unsigned, double* out) {
  const auto& choice = *static_cast<const ProbitChoice*>(data);
  double p = std::exp(-0.5 * x[0] * x[0]) * M_1_SQRT_2PI;
  for (int j = 0; j < choice.n; ++j) {
    if (j == choice.winner) continue;
    const double z = x[0] + (choice.utility[choice.winner] - choice.utility[j]);
    p *= 0.5 * std::erfc(-z * M_SQRT1_2);
  }
  out[0] = p;
  return 0;
}

}

double multinomial_probit_probability(const double* utility, int n, int winner) {
  if (n == 1) return 1.0;
  if (n == 2)
    return 0.5 * std::erfc(-(utility[winner] - utility[1 - winner]) / 2.0);

  // P(winner) = integral phi(x) product_{j != winner} Phi(x + u_winner - u_j) dx.
  // The omitted standard-normal tails have total mass below 4e-33.
  ProbitChoice choice{utility, n, winner};
  const double lower = -12.0, upper = 12.0;
  double p, error;
  hcubature(probit_integrand, &choice, 1, &lower, &upper, 4095,
            1e-24, 1e-10, &p, &error);
  return p;
}

// [[Rcpp::export]]
Rcpp::NumericVector c_multinomial_probit_probs(Rcpp::NumericVector utility) {
  const int n = utility.size();
  Rcpp::NumericVector probs(n);
  double total = 0.0;
  for (int i = 0; i < n; ++i) {
    probs[i] = multinomial_probit_probability(utility.begin(), n, i);
    total += probs[i];
  }
  return probs / total;
}
