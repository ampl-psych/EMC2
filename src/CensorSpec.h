#ifndef CENSORSPEC_H
#define CENSORSPEC_H

#include <vector>

struct CensorSpec {
  std::vector<int>    idx_L;   // lower-censored rows
  std::vector<int>    idx_U;   // upper-censored rows
  std::vector<int>    idx_B;   // both-censored rows
  std::vector<double> LC;      // lower bounds (length n_trials, or empty)
  std::vector<double> UC;      // upper bounds (length n_trials, or empty)

  bool any() const { return !idx_L.empty() || !idx_U.empty() || !idx_B.empty(); }
};

CensorSpec make_censor_spec(const Rcpp::DataFrame& data,
                            int n_trials);

#endif //CENSORSPEC_H
