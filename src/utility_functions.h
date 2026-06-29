#ifndef utility_h
#define utility_h

#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <functional>
#include "ParamTable.h"
using namespace Rcpp;

// ---------------------------------------------------------------------------
// -ffast-math-safe finite/Inf/NaN checks. Under -ffast-math the compiler may
// constant-fold std::isinf/std::isfinite and never inspect the bit pattern, so
// for runtime VALUES we use R's bit-pattern macros (R_FINITE, ISNAN) which the
// compiler cannot optimise away.
//   emc2_isfinite(x)  ->  R_FINITE(x)            (0 for Inf AND NaN)
//   emc2_isinf(x)     ->  !R_FINITE(x) && !ISNAN(x)
//   emc2_isnan(x)     ->  ISNAN(x)
// ---------------------------------------------------------------------------
static inline bool emc2_isfinite(double x) { return (bool)R_FINITE(x); }
static inline bool emc2_isinf(double x)    { return !R_FINITE(x) && !ISNAN(x); }
static inline bool emc2_isnan(double x)    { return (bool)ISNAN(x); }


// Custom hash for a vector<double>
struct RowHash {
  std::size_t operator()(const std::vector<double> &v) const {
    std::size_t seed = 0;
    std::hash<double> hash_double;
    for (double d : v) {
      // A standard hash combination approach (based on boost::hash_combine)
      seed ^= hash_double(d) + 0x9e3779b97f4a7c16ULL + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// Custom equality for a vector<double>
struct RowEqual {
  bool operator()(const std::vector<double> &a, const std::vector<double> &b) const {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
      if (a[i] != b[i]) return false;
    }
    return true;
  }
};

LogicalVector contains(CharacterVector sv, std::string txt);

NumericVector vector_pow(NumericVector x1, NumericVector x2);

NumericVector pnorm_multiple(NumericVector x);

LogicalVector contains_multiple(CharacterVector sv, CharacterVector inputs);

NumericMatrix submat_rcpp_col(NumericMatrix X, LogicalVector condition);

NumericMatrix submat_rcpp_col_by_names(NumericMatrix X, CharacterVector cols);

NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition);

NumericVector c_expand(NumericVector x1, IntegerVector expand);

LogicalVector c_bool_expand(LogicalVector x1, IntegerVector expand);

NumericVector c_add_vectors(NumericVector x1, NumericVector x2);

CharacterVector c_add_charvectors(CharacterVector x, CharacterVector y);

LogicalVector duplicated_matrix(Rcpp::NumericMatrix x);

IntegerVector cumsum_logical(LogicalVector x);

IntegerVector which_rcpp(LogicalVector x);

void lr_all(std::vector<int>& ok, int n_side);

#endif

