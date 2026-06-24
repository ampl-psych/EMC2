#include <Rcpp.h>
#include <algorithm>
#include "CensorSpec.h"

using namespace Rcpp;

CensorSpec make_censor_spec(const DataFrame& data, int n_trials)
{
  CensorSpec censor;

  std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(data.names());
  const bool has_LC_col      = std::find(names.begin(), names.end(), "LC")          != names.end();
  const bool has_UC_col      = std::find(names.begin(), names.end(), "UC")          != names.end();
  const bool has_missingness = std::find(names.begin(), names.end(), "missingness") != names.end();

  if (!has_missingness) return censor;

  // Raw pointers — kept alive by the Rcpp wrappers below
  IntegerVector missingness_tmp = data["missingness"];
  const int* missingness_ptr    = INTEGER(missingness_tmp);

  NumericVector lc_tmp, uc_tmp;
  if (has_LC_col) {
    lc_tmp = data["LC"];
    censor.LC.assign(REAL(lc_tmp), REAL(lc_tmp) + n_trials);
    censor.idx_L.reserve(n_trials);
  }
  if (has_UC_col) {
    uc_tmp = data["UC"];
    censor.UC.assign(REAL(uc_tmp), REAL(uc_tmp) + n_trials);
    censor.idx_U.reserve(n_trials);
    censor.idx_B.reserve(n_trials);
  }

  for (int i = 0; i < n_trials; ++i) {
    switch (missingness_ptr[i]) {
    case 1: censor.idx_L.push_back(i); break;
    case 2: censor.idx_U.push_back(i); break;
    case 3: censor.idx_B.push_back(i); break;
    default: break;
    }
  }

  // Validate consistency
  if (!censor.idx_L.empty() && censor.LC.empty())
    Rcpp::stop("missingness contains lower-censored rows (1 or 3) but no LC column found in data");
  if (!censor.idx_U.empty() && censor.UC.empty())
    Rcpp::stop("missingness contains upper-censored rows (2 or 3) but no UC column found in data");

  return censor;
}
