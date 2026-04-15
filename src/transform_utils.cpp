// transform_utils.cpp
#include "transform_utils.h"

using namespace Rcpp;
using std::string;

// ---------------------
// make_transform_specs (NumericMatrix)
// ---------------------

std::vector<TransformSpec> make_transform_specs_matrix(const NumericMatrix& pars, const List& transform) {
    CharacterVector func_charvec = transform["func"];
    NumericVector   lower_numvec = transform["lower"];
    NumericVector   upper_numvec = transform["upper"];

    // name -> code
    std::unordered_map<string,TransformCode> codeMap;
    codeMap.reserve(func_charvec.size());
    {
      CharacterVector fnames = func_charvec.names();
      for (int i = 0; i < func_charvec.size(); i++) {
        string name = as<string>(fnames[i]);
        string f    = as<string>(func_charvec[i]);
        TransformCode code = IDENTITY;
        if (f == "exp")      code = EXP;
        else if (f == "pnorm") code = PNORM;
        codeMap.emplace(name, code);
      }
    }

    // name -> (lower, upper)
    std::unordered_map<string,std::pair<double,double>> boundMap;
    {
      CharacterVector ln = lower_numvec.names();
      for (int i = 0; i < lower_numvec.size(); i++) {
        string nm = as<string>(ln[i]);
        boundMap[nm].first = lower_numvec[i];
      }
      CharacterVector un = upper_numvec.names();
      for (int i = 0; i < upper_numvec.size(); i++) {
        string nm = as<string>(un[i]);
        boundMap[nm].second = upper_numvec[i];
      }
    }

    int ncol = pars.ncol();
    std::vector<TransformSpec> specs(ncol);
    CharacterVector cparnames = colnames(pars);

    for (int j = 0; j < ncol; j++) {
      string colname = as<string>(cparnames[j]);
      TransformSpec sp;
      sp.col_idx = j;

      auto itc = codeMap.find(colname);
      sp.code = (itc != codeMap.end()) ? itc->second : IDENTITY;

      auto itb = boundMap.find(colname);
      if (itb != boundMap.end()) {
        sp.lower = itb->second.first;
        sp.upper = itb->second.second;
      } else {
        sp.lower = 0.0;
        sp.upper = 1.0;
      }
      specs[j] = sp;
    }

    return specs;
  }


// ---------------------
// make_transform_specs_pt
// ---------------------

std::vector<TransformSpec> make_transform_specs_pt(const ParamTable& pt, const List& transform) {
    CharacterVector func_charvec = transform["func"];
    NumericVector lower_numvec   = transform["lower"];
    NumericVector upper_numvec   = transform["upper"];

    std::unordered_map<string, TransformCode> codeMap;
    codeMap.reserve(func_charvec.size());
    {
      CharacterVector fnames = func_charvec.names();
      for (int i = 0; i < func_charvec.size(); ++i) {
        string name = as<string>(fnames[i]);
        string f    = as<string>(func_charvec[i]);
        TransformCode code = IDENTITY;
        if (f == "exp")      code = EXP;
        else if (f == "pnorm") code = PNORM;
        codeMap.emplace(name, code);
      }
    }

    std::unordered_map<string, std::pair<double,double>> boundMap;
    {
      CharacterVector ln = lower_numvec.names();
      for (int i = 0; i < lower_numvec.size(); ++i) {
        boundMap[ as<string>(ln[i]) ].first = lower_numvec[i];
      }
      CharacterVector un = upper_numvec.names();
      for (int i = 0; i < upper_numvec.size(); ++i) {
        boundMap[ as<string>(un[i]) ].second = upper_numvec[i];
      }
    }

    const int n_active = pt.n_params();
    std::vector<TransformSpec> specs;
    specs.reserve(n_active);

    for (int k = 0; k < n_active; ++k) {
      int base_idx = pt.active_cols[k];
      string colname = as<string>(pt.base_names[base_idx]);

      TransformSpec sp;
      sp.col_idx = base_idx;

      auto itc = codeMap.find(colname);
      sp.code = (itc != codeMap.end()) ? itc->second : IDENTITY;

      auto itb = boundMap.find(colname);
      if (itb != boundMap.end()) {
        sp.lower = itb->second.first;
        sp.upper = itb->second.second;
      } else {
        sp.lower = 0.0;
        sp.upper = 1.0;
      }

      specs.push_back(sp);
    }

    return specs;
  }


// ---------------------
// c_do_transform (NumericMatrix)
// ---------------------

NumericMatrix c_do_transform_matrix(NumericMatrix pars, const std::vector<TransformSpec>& specs) {
  int nrow = pars.nrow();

  for (size_t j = 0; j < specs.size(); j++) {
    const TransformSpec& sp = specs[j];
    int          col_idx = sp.col_idx;
    TransformCode c      = sp.code;
    double        lw     = sp.lower;
    double        up     = sp.upper;

    switch (c) {
    case EXP: {
      for (int i = 0; i < nrow; i++) {
      pars(i, col_idx) = lw + std::exp(pars(i, col_idx));
    }
      break;
    }
    case PNORM: {
      double range = up - lw;
      for (int i = 0; i < nrow; i++) {
        pars(i, col_idx) = lw +
          range * R::pnorm(pars(i, col_idx), 0.0, 1.0, 1, 0);
      }
      break;
    }
    case IDENTITY:
    default:
      break;
    }
  }
  return pars;
}


Rcpp::LogicalVector c_do_bound_pt(const ParamTable& pt,
                                  const std::vector<BoundSpec>& specs)
{
  using Rcpp::LogicalVector;
  using Rcpp::NumericMatrix;

  const NumericMatrix& base = pt.base;
  const int nrows = base.nrow();

  LogicalVector result(nrows, true);

  for (std::size_t j = 0; j < specs.size(); ++j) {
    const BoundSpec& bs = specs[j];
    const int col_idx   = bs.col_idx;   // MUST be a base-column index
    const double min_v  = bs.min_val;
    const double max_v  = bs.max_val;
    const bool has_exc  = bs.has_exception;
    const double exc_val= bs.exception_val;

    for (int i = 0; i < nrows; ++i) {
      const double val = base(i, col_idx);
      bool ok = (val > min_v && val < max_v);
      if (!ok && has_exc) {
        ok = (val == exc_val);
      }
      if (result[i] && !ok) {
        result[i] = false;
      }
    }
  }

  return result;
}

// ---------------------
// c_do_transform_pt (ParamTable version)
// ---------------------

void c_do_transform_pt(ParamTable& pt,
                  const std::vector<TransformSpec>& specs)
{
  const int nrow = pt.n_trials;

  for (size_t j = 0; j < specs.size(); ++j) {
    const TransformSpec& sp = specs[j];
    const int col_idx       = sp.col_idx;
    const TransformCode c   = sp.code;
    const double lw         = sp.lower;
    const double up         = sp.upper;

    double* col = &pt.base(0, col_idx);

    switch (c) {
    case EXP: {
      for (int i = 0; i < nrow; ++i) {
      col[i] = lw + std::exp(col[i]);
    }
      break;
    }
    case PNORM: {
      const double range = up - lw;
      for (int i = 0; i < nrow; ++i) {
        col[i] = lw + range *
          R::pnorm(col[i], 0.0, 1.0, 1, 0);
      }
      break;
    }
    case IDENTITY:
    default:
      break;
    }
  }
}

// Filters for trend parameter names corresponding to trends that need to be applied pre-map
std::vector<TransformSpec> filter_specs_by_param_set(
    const ParamTable& pt,
    const std::vector<TransformSpec>& full_specs,
    const std::unordered_set<std::string>& allowed)
  {
  std::vector<TransformSpec> out;
  out.reserve(full_specs.size());

  for (const auto& sp : full_specs) {
    int base_idx = sp.col_idx;
    std::string nm = Rcpp::as<std::string>(pt.base_names[base_idx]);
    if (allowed.find(nm) != allowed.end()) {
      out.push_back(sp);
    }
  }
  return out;
}

// Filters for all other parameter names
std::vector<TransformSpec> complement_specs_for_premap(
    const ParamTable& pt,
    const std::vector<TransformSpec>& full_specs,
    const std::unordered_set<std::string>& premap_set)
{
  std::vector<TransformSpec> out;
  out.reserve(full_specs.size());

  for (const auto& sp : full_specs) {
    int base_idx = sp.col_idx;
    std::string nm = Rcpp::as<std::string>(pt.base_names[base_idx]);
    if (premap_set.find(nm) == premap_set.end()) {
      out.push_back(sp);
    }
  }
  return out;
}


// Same logic as above, but with a ParamTable instead of NumericMatrix
std::vector<BoundSpec> make_bound_specs_pt(Rcpp::NumericMatrix minmax,
                                           Rcpp::CharacterVector minmax_colnames,
                                           const ParamTable& pt,
                                           Rcpp::List bound)
{
  using namespace Rcpp;
  using std::string;

  // 1) Build a map from param-name -> base column index in ParamTable
  //    (or just use pt.base_index_for(...))
  // Here we go directly via pt.base_index_for for simplicity.

  // 2) Build a map from param-name -> exception value
  bool has_exception =
    bound.containsElementNamed("exception") &&
    !Rf_isNull(bound["exception"]);

    std::unordered_map<string, double> exceptionMap;
    if (has_exception) {
      NumericVector except_vec = bound["exception"];
      CharacterVector except_names = except_vec.names();
      for (int i = 0; i < except_vec.size(); ++i) {
        string nm = as<string>(except_names[i]);
        exceptionMap[nm] = except_vec[i];
      }
    }

    // 3) Create BoundSpec for each column in minmax
    const int ncols = minmax_colnames.size();
    std::vector<BoundSpec> specs(ncols);

    for (int j = 0; j < ncols; ++j) {
      string var_name = as<string>(minmax_colnames[j]);

      // find base-column index for this parameter
      int base_idx = pt.base_index_for(var_name);  // will Rcpp::stop if unknown

      BoundSpec s;
      s.col_idx = base_idx;           // index into pt.base
      s.min_val = minmax(0, j);
      s.max_val = minmax(1, j);

      auto it = exceptionMap.find(var_name);
      if (it != exceptionMap.end()) {
        s.has_exception = true;
        s.exception_val = it->second;
      } else {
        s.has_exception = false;
        s.exception_val = NA_REAL;  // same as before
      }

      specs[j] = s;
    }

    return specs;
}
