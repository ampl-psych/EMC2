#pragma once

#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
using Rcpp::_;

// View-based ParamTable: one base matrix + active column indices
// Invariants:
// - base contains all parameters (never shrinks)
// - active_cols is a view into base (ordering may change)
// - name_to_base_idx maps ALL base columns, active or not
struct ParamTable {
  Rcpp::NumericMatrix base;     // underlying matrix (all columns)
  Rcpp::CharacterVector base_names; // colnames for base
  int n_trials = 0;

  // I was thinking of including some logic that keeps track of which columns are still 'active'
  // but that's not actually used at this point..
  std::vector<int> active_cols;  // indices into base/base_names
  std::vector<char> is_active;  // size = base.ncol()

  // fast look-up map
  std::unordered_map<std::string,int> name_to_base_idx;

  // keep track of parameters that have an intercept-only design
  std::vector<char> design_is_self_intercept;

  ParamTable() = default;

  ParamTable(Rcpp::NumericMatrix base_,
             Rcpp::CharacterVector names_)
    : base(base_), base_names(names_) {
    n_trials = base_.nrow();
    const int p = base_.ncol();
    active_cols.resize(p);
    std::iota(active_cols.begin(), active_cols.end(), 0);
    Rcpp::colnames(base) = base_names;
    rebuild_name_map();
  }

  void rebuild_name_map() {
    name_to_base_idx.clear();
    const int p = base_names.size();
    name_to_base_idx.reserve(p);
    for (int j = 0; j < p; ++j) {
      name_to_base_idx[ Rcpp::as<std::string>(base_names[j]) ] = j;
    }
  }

  int base_index_for(const std::string& nm) const {
    auto it = name_to_base_idx.find(nm);
    if (it == name_to_base_idx.end()) {
      Rcpp::stop("ParamTable: Unknown parameter '%s'", nm.c_str());
    }
    return it->second;
  }
  int n_params() const { return (int)active_cols.size(); }

  // Name of j-th active column
  Rcpp::String name_at(int j) const {
    return base_names[ active_cols[j] ];
  }

  // Column view by active index (no deep copy)
  Rcpp::NumericVector column_by_active_index(int j) const {
    return base(_, active_cols[j]);  // view
  }

  bool is_active_base_idx(int base_idx) const {
    return std::find(active_cols.begin(), active_cols.end(), base_idx) != active_cols.end();
  }

  // Column view by parameter name
  Rcpp::NumericVector column_by_name(const std::string& nm) const {
    int base_idx = base_index_for(nm); // O(1) from map
    // Ensure it's currently active; if you need this check:
    // (If you know all base cols are always active, you can skip this loop.)
    // for (int j = 0; j < (int)active_cols.size(); ++j) {
      // if (active_cols[j] == base_idx) {
        return base(_, base_idx);  // view
      // }
    // }
    // Rcpp::stop("ParamTable: parameter '%s' not active", nm.c_str());
  }

  // Assign into existing column (writes into base)
  void set_column_by_name(const std::string& nm,
                          const Rcpp::NumericVector& col) {
    int base_idx = base_index_for(nm);
    // Optional: ensure active
    // for (int j = 0; j < (int)active_cols.size(); ++j) {
      // if (active_cols[j] == base_idx) {
        base(_, base_idx) = col;
        return;
      // }
    // }
    // Rcpp::stop("ParamTable: parameter '%s' not active", nm.c_str());
  }

  // Drop columns by name: only adjust active_cols
  // void drop(const Rcpp::CharacterVector& drop_names) {
  //   std::unordered_set<std::string> drop_set;
  //   for (int i = 0; i < drop_names.size(); ++i)
  //     drop_set.insert(Rcpp::as<std::string>(drop_names[i]));
  //
  //   std::vector<int> new_active;
  //   new_active.reserve(active_cols.size());
  //   for (int idx : active_cols) {
  //     std::string nm = Rcpp::as<std::string>(base_names[idx]);
  //     if (!drop_set.count(nm)) new_active.push_back(idx);
  //   }
  //   active_cols.swap(new_active);
  // }

  Rcpp::NumericMatrix materialize_by_param_names(const Rcpp::CharacterVector& param_names) const {
    using namespace Rcpp;
    const int k = param_names.size();
    NumericMatrix out(n_trials, k);
    CharacterVector out_names(k);

    for (int j = 0; j < k; ++j) {
      std::string nm = as<std::string>(param_names[j]);
      int base_idx = base_index_for(nm);  // throws if unknown

      // Ensure this column is active (if you care about active_cols)
      bool is_active = false;
      for (int a : active_cols) {
        if (a == base_idx) {
          is_active = true;
          break;
        }
      }
      if (!is_active) {
        stop("ParamTable::materialize_by_param_names: parameter '%s' not active", nm.c_str());
      }

      // Copy base(:, base_idx) → out(:, j)
      double* out_col        = &out(0, j);
      const double* base_col = &base(0, base_idx);
      for (int r = 0; r < n_trials; ++r) {
        out_col[r] = base_col[r];
      }
      out_names[j] = base_names[base_idx];
    }

    colnames(out) = out_names;
    return out;
  }

  // Materialize matrix (n_trials x n_active) with all columns
  Rcpp::NumericMatrix materialize() const {
    const int p = (int)active_cols.size();
    Rcpp::NumericMatrix out(n_trials, p);
    Rcpp::CharacterVector out_names(p);

    for (int j = 0; j < p; ++j) {
      int base_j = active_cols[j];
      // avoid Rcpp Sugar copies like out(_, j) = base(_, base_j);
      double* out_col = &out(0, j);
      const double* base_col = &base(0, base_j);
      for (int r = 0; r < n_trials; ++r) {
        out_col[r] = base_col[r];
      }
      out_names[j] = base_names[base_j];
    }
    Rcpp::colnames(out) = out_names;
    return out;
  }

  static ParamTable from_p_types(int n_trials,
                                 const Rcpp::CharacterVector& p_types) {
    const int p = p_types.size();
    Rcpp::NumericMatrix base(n_trials, p);
    base.fill(0.0);
    Rcpp::CharacterVector names = Rcpp::clone(p_types);
    return ParamTable(base, names);
  }

  static ParamTable from_p_vector_and_designs(const Rcpp::NumericVector& p_vector,
                                              const Rcpp::List& designs,
                                              int n_trials) {
    using std::string;
    using Rcpp::as;
    using Rcpp::CharacterVector;
    using Rcpp::NumericMatrix;
    using Rcpp::NumericVector;

    // 1) Collect all unique parameter names in a stable order:
    //    (a) names(p_vector)
    //    (b) names(designs)
    //    (c) colnames(designs[[i]])
    std::vector<string> names_vec;
    names_vec.reserve(p_vector.size() + designs.size() * 2);
    std::unordered_set<string> seen;

    // (a) names(p_vector)
    CharacterVector pv_names = p_vector.names();
    for (int i = 0; i < pv_names.size(); ++i) {
      string nm = as<string>(pv_names[i]);
      if (!seen.count(nm)) {
        seen.insert(nm);
        names_vec.push_back(nm);
      }
    }

    // (b) names(designs)
    CharacterVector design_names = designs.names();
    for (int i = 0; i < design_names.size(); ++i) {
      string nm = as<string>(design_names[i]);
      if (!seen.count(nm)) {
        seen.insert(nm);
        names_vec.push_back(nm);
      }
    }

    // (c) colnames of each design matrix
    for (int i = 0; i < designs.size(); ++i) {
      if (designs[i] == R_NilValue) continue;
      NumericMatrix dm = designs[i];
      CharacterVector cn = Rcpp::colnames(dm);
      for (int j = 0; j < cn.size(); ++j) {
        string nm = as<string>(cn[j]);
        if (!seen.count(nm)) {
          seen.insert(nm);
          names_vec.push_back(nm);
        }
      }
    }

    // 2) Build base matrix and names vector
    const int p = static_cast<int>(names_vec.size());
    NumericMatrix base(n_trials, p);
    base.fill(0.0);

    CharacterVector base_names(p);
    for (int j = 0; j < p; ++j) {
      base_names[j] = names_vec[j];
    }
    Rcpp::colnames(base) = base_names;

    // 3) Build map name -> value from p_vector
    std::unordered_map<string, double> pval;
    pval.reserve(p_vector.size());
    for (int i = 0; i < p_vector.size(); ++i) {
      string nm = as<string>(pv_names[i]);
      pval[nm] = p_vector[i];
    }

    // 4) Fill columns that appear in p_vector with that scalar value, others remain 0
    for (int j = 0; j < p; ++j) {
      string nm = names_vec[j];
      auto it = pval.find(nm);
      if (it != pval.end()) {
        double val = it->second;
        for (int r = 0; r < n_trials; ++r) {
          base(r, j) = val;
        }
      }
      // else: not in p_vector -> keep zeros
    }

    // 5) Construct ParamTable
    // return ParamTable(base, base_names);

    // 5) Construct ParamTable
    ParamTable pt(base, base_names);
    pt.n_trials = n_trials;  // (already set in ctor, but fine to be explicit)

    // cache "self-intercept-only" designs
    const int n_designs = designs.size();
    pt.design_is_self_intercept.assign(n_designs, 0);

    for (int i = 0; i < n_designs; ++i) {
      if (designs[i] == R_NilValue) continue;

      NumericMatrix dm = designs[i];

      // Must have the right number of rows
      if (dm.nrow() != n_trials) continue;

      // Exactly 1 column
      if (dm.ncol() != 1) continue;

      CharacterVector cn = Rcpp::colnames(dm);
      if (cn.size() != 1) continue;

      // Column name must match the design/output name
      string out_name  = as<string>(design_names[i]);
      string coef_name = as<string>(cn[0]);
      if (coef_name != out_name) continue;

      // Column must be all ones
      bool all_ones = true;
      const double* dcol = &dm(0, 0);
      for (int r = 0; r < n_trials; ++r) {
        if (dcol[r] != 1.0) {
          all_ones = false;
          break;
        }
      }

      if (all_ones) {
        pt.design_is_self_intercept[i] = 1;
      }
    }

    return pt;
  }

  // Map parameters from designs into this table, ignoring trends.
  // Uses the *current* trial-by-trial values in ParamTable::base as
  // coefficients for design columns.
  //
  // designs: named list of design matrices (possibly compressed).
  //          Each matrix may have an "expand" attribute:
  //          - if present and length == n_trials: use it to expand rows
  //          - else, require nrow(design) == n_trials.
  // include_param: optional logical mask; if length == length(designs),
  //                decides which entries to map; if empty, all TRUE.
  void map_from_designs(const Rcpp::List& designs,
                        const Rcpp::LogicalVector& include_param = Rcpp::LogicalVector()) {
    using namespace Rcpp;
    using std::string;

    const int n_params = designs.size();
    if (n_params == 0) return;

    CharacterVector design_names = designs.names();
    if (design_names.size() != n_params) {
      stop("ParamTable::map_from_designs: designs must be a named list");
    }

    LogicalVector use =
      (include_param.size() == n_params)
      ? include_param
    : LogicalVector(n_params, true);

    const int T = n_trials;

    for (int i = 0; i < n_params; ++i) {
      if (!use[i]) continue;
      if (designs[i] == R_NilValue) continue;

      // skip self-intercept-only designs -- multiplying by 1 is pointless
      if (!design_is_self_intercept.empty() &&
          i < (int)design_is_self_intercept.size() &&
          design_is_self_intercept[i]) {
        // Rprintf("Skipping column as this is a self-intercept \n");
        continue;
      }

      const string out_name = as<string>(design_names[i]);

      // --- Resolve output column ---
      auto out_it = name_to_base_idx.find(out_name);
      if (out_it == name_to_base_idx.end()) {
        stop("ParamTable::map_from_designs: output parameter '%s' not in table",
             out_name.c_str());
      }
      const int out_idx = out_it->second;

      NumericMatrix design = designs[i];
      if (design.nrow() != T) {
        stop("ParamTable::map_from_designs: design for '%s' must have n_trials rows",
             out_name.c_str());
      }

      const int K = design.ncol();
      CharacterVector coef_names = colnames(design);

      // --- Resolve coefficient indices ---
      std::vector<int> coef_idx(K, -1);
      bool uses_self = false;

      for (int j = 0; j < K; ++j) {
        auto it = name_to_base_idx.find(as<string>(coef_names[j]));
        if (it == name_to_base_idx.end()) continue;
        coef_idx[j] = it->second;
        if (coef_idx[j] == out_idx) uses_self = true;
      }

      // --- Preserve self column if needed ---
      std::vector<double> self_copy;
      if (uses_self) {
        self_copy.resize(T);
        const double* src = &base(0, out_idx);
        std::copy(src, src + T, self_copy.begin());
      }

      // --- Clear output ---
      double* out = &base(0, out_idx);
      std::fill(out, out + T, 0.0);

      // --- Main accumulation loop ---
      for (int j = 0; j < K; ++j) {
        int cidx = coef_idx[j];
        if (cidx < 0) continue;

        const double* coef =
          (cidx == out_idx && uses_self)
          ? self_copy.data()
            : &base(0, cidx);

        const double* d = &design(0, j);

        for (int r = 0; r < T; ++r) {
          double v = coef[r] * d[r];
          // if (!std::isfinite(v)) continue;
          out[r] += v;
        }
      }
    }
  }

  // Zero the entire base matrix
  void reset_base_to_zero() {
    const int n = n_trials;
    const int p = base.ncol();
    for (int j = 0; j < p; ++j) {
      double* col = &base(0, j);
      for (int r = 0; r < n; ++r) {
        col[r] = 0.0;
      }
    }
  }

  // Fill columns corresponding to names(p_vector) with that scalar value
  void fill_from_p_vector(const Rcpp::NumericVector& p_vector) {
    using std::string;
    using Rcpp::as;

    Rcpp::CharacterVector pv_names = p_vector.names();
    const int n_names = pv_names.size();

    for (int i = 0; i < n_names; ++i) {
      string nm = as<string>(pv_names[i]);
      auto it = name_to_base_idx.find(nm);
      if (it == name_to_base_idx.end()) {
        // p_vector may contain names that weren't used in designs, ignore them
        continue;
      }
      int j = it->second;
      double val = p_vector[i];
      double* col = &base(0, j);
      for (int r = 0; r < n_trials; ++r) {
        col[r] = val;
      }
    }
  }
};


// Helper: all parameter names in pt that are NOT in premap_trend_params
std::unordered_set<std::string>
  make_non_premap_param_set(const ParamTable& pt,
                            const std::unordered_set<std::string>& premap_trend_params);

// Helper: all parameter names in pt that are NOT in premap_trend_params OR in pretransform_trend_params
std::unordered_set<std::string>
  make_non_premap_pretransform_param_set(const ParamTable& pt,
                                         const std::unordered_set<std::string>& premap_trend_params,
                                         const std::unordered_set<std::string>& pretransform_trend_params);
