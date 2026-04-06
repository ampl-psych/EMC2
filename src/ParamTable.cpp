#include "ParamTable.h"


Rcpp::CharacterVector names_excluding(const Rcpp::CharacterVector& names,
                                     std::initializer_list<const std::unordered_set<std::string>*> excludes)
{
  using namespace Rcpp;
  using std::string;

  std::vector<string> keep;
  keep.reserve(names.size());

  for (int i = 0; i < names.size(); ++i) {
    string nm = as<string>(names[i]);

    bool skip = false;
    for (auto set_ptr : excludes) {
      if (!set_ptr) continue;
      const auto& s = *set_ptr;
      if (s.find(nm) != s.end()) {
        skip = true;
        break;
      }
    }
    if (!skip) {
      keep.push_back(std::move(nm));
    }
  }

  CharacterVector out(keep.size());
  for (int i = 0; i < (int)keep.size(); ++i) {
    out[i] = keep[i];
  }
  return out;
}

// All parameter names in pt that are NOT in any of the exclude sets.
//
// Usage examples:
//   auto all_but_premap       = param_names_excluding(pt, { &premap_set });
//   auto all_but_premap_pretr = param_names_excluding(pt, { &premap_set, &pretransform_set });
//
std::unordered_set<std::string> param_names_excluding(const ParamTable& pt,
                                                             std::initializer_list<const std::unordered_set<std::string>*> excludes)
{
  using namespace Rcpp;
  using std::string;

  // reuse the generic logic but on pt.base_names
  CharacterVector base_names = pt.base_names;
  CharacterVector kept = names_excluding(base_names, excludes);

  std::unordered_set<string> out;
  out.reserve(kept.size());
  for (int i = 0; i < kept.size(); ++i) {
    out.insert(as<string>(kept[i]));
  }
  return out;
}


// Little helper for constants
Rcpp::NumericMatrix add_constants_columns(Rcpp::NumericMatrix p_matrix,
                                          Rcpp::NumericVector constants) {
  const int n_rows = p_matrix.nrow();
  const int p_old  = p_matrix.ncol();
  const int p_add  = constants.size();
  const int p_new  = p_old + p_add;

  // allocate new matrix
  Rcpp::NumericMatrix out(n_rows, p_new);

  // copy existing data
  for (int j = 0; j < p_old; ++j) {
    for (int i = 0; i < n_rows; ++i) {
      out(i, j) = p_matrix(i, j);
    }
  }

  // fill new columns with constants (same value for all rows)
  for (int k = 0; k < p_add; ++k) {
    double val = constants[k];
    int j = p_old + k;
    for (int i = 0; i < n_rows; ++i) {
      out(i, j) = val;
    }
  }

  // set column names: old names + names(constants)
  Rcpp::CharacterVector old_names = colnames(p_matrix);
  Rcpp::CharacterVector const_names = constants.names();
  Rcpp::CharacterVector new_names(p_new);

  for (int j = 0; j < p_old; ++j) {
    new_names[j] = old_names[j];
  }
  for (int k = 0; k < p_add; ++k) {
    new_names[p_old + k] = const_names[k];
  }

  colnames(out) = new_names;

  return out;
}
