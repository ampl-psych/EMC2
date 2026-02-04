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
    return ParamTable(base, base_names);
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
  // Map parameters from designs into this table, assuming *every*
  // non-NULL design matrix has an "expand" attribute of length n_trials.
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
  // void map_from_designs(const Rcpp::List& designs,
  //                       const Rcpp::LogicalVector& include_param = Rcpp::LogicalVector()) {
  //   using namespace Rcpp;
  //   using std::string;
  //
  //   const int n_params = designs.size();
  //   if (n_params == 0) return;
  //
  //   CharacterVector design_names = designs.names();
  //   if (design_names.size() != n_params) {
  //     stop("ParamTable::map_from_designs: designs must be a named list");
  //   }
  //
  //   // Decide which parameters to map
  //   LogicalVector use;
  //   if (include_param.size() == n_params) {
  //     use = include_param;
  //   } else {
  //     use = LogicalVector(n_params, true); // default: map all
  //   }
  //
  //   const int n_trials_local = n_trials;
  //
  //   for (int i = 0; i < n_params; ++i) {
  //     if (!use[i]) continue;
  //
  //     string out_pname = as<string>(design_names[i]);
  //
  //     if (designs[i] == R_NilValue) {
  //       // No design for this parameter → leave as-is
  //       continue;
  //     }
  //
  //     NumericMatrix cur_design = designs[i];
  //     const int n_rows_comp = cur_design.nrow();
  //     CharacterVector cur_names = colnames(cur_design);
  //     const int n_cols = cur_design.ncol();
  //
  //     // Output column index (must exist)
  //     auto out_it = name_to_base_idx.find(out_pname);
  //     if (out_it == name_to_base_idx.end()) {
  //       stop("ParamTable::map_from_designs: output parameter '%s' not in table",
  //            out_pname.c_str());
  //     }
  //     const int out_idx = out_it->second;
  //
  //     // --- Handle "expand" attribute: may be NULL ---
  //     IntegerVector expand_idx;
  //     {
  //       SEXP exp_attr = cur_design.attr("expand");
  //       if (exp_attr == R_NilValue) {
  //         // No expand attribute: require full-length design
  //         if (n_rows_comp != n_trials_local) {
  //           stop("ParamTable::map_from_designs: design for '%s' has no 'expand' "
  //                  "attribute and nrow(design) != n_trials (%d != %d)",
  //                  out_pname.c_str(), n_rows_comp, n_trials_local);
  //         }
  //         // identity expand: trial r uses row r
  //         expand_idx = IntegerVector(n_trials_local);
  //         for (int r = 0; r < n_trials_local; ++r) {
  //           expand_idx[r] = r + 1;  // R is 1-based
  //         }
  //       } else {
  //         expand_idx = exp_attr;   // will throw if not integer-like
  //         if (expand_idx.size() != n_trials_local) {
  //           stop("ParamTable::map_from_designs: length(expand) != n_trials for parameter '%s'",
  //                out_pname.c_str());
  //         }
  //       }
  //     }
  //     // --- Pre-compute coefficient indices and detect self-coefficient use ---
  //     std::vector<int> coef_indices(n_cols, -1);
  //     bool uses_self = false;
  //     for (int j = 0; j < n_cols; ++j) {
  //       string col_name = as<string>(cur_names[j]);
  //       auto coef_it = name_to_base_idx.find(col_name);
  //       if (coef_it == name_to_base_idx.end()) {
  //         coef_indices[j] = -1;  // "missing" => contributes 0
  //         continue;
  //       }
  //       int coef_idx = coef_it->second;
  //       coef_indices[j] = coef_idx;
  //       if (coef_idx == out_idx) {
  //         uses_self = true;
  //       }
  //     }
  //
  //     // If the design uses out_pname itself as a coefficient, copy its
  //     // original values BEFORE we clear the output column.
  //     std::vector<double> self_coef;
  //     if (uses_self) {
  //       self_coef.resize(n_trials_local);
  //       const double* src = &base(0, out_idx);
  //       for (int r = 0; r < n_trials_local; ++r) {
  //         self_coef[r] = src[r];
  //       }
  //     }
  //
  //     // Clear output column and get raw pointer to it
  //     double* out_ptr = &base(0, out_idx);
  //     for (int r = 0; r < n_trials_local; ++r) {
  //       out_ptr[r] = 0.0;
  //     }
  //
  //     // Reusable buffer for expanded design column
  //     NumericVector design_col(n_trials_local);
  //
  //     // --- Main loop over design columns ---
  //     for (int j = 0; j < n_cols; ++j) {
  //       int coef_idx = coef_indices[j];
  //       if (coef_idx < 0) continue;  // no such coefficient -> contributes 0
  //
  //       // Choose coefficient pointer:
  //       const double* coef_ptr = nullptr;
  //       if (coef_idx == out_idx && uses_self) {
  //         // Use the saved original column if it's the self-coefficient
  //         coef_ptr = self_coef.data();
  //       } else {
  //         coef_ptr = &base(0, coef_idx);
  //       }
  //
  //       // Expand compressed design column j into design_col
  //       for (int r = 0; r < n_trials_local; ++r) {
  //         int row_idx = expand_idx[r] - 1; // R 1-based → C++ 0-based
  //         if (row_idx < 0 || row_idx >= n_rows_comp) {
  //           stop("ParamTable::map_from_designs: 'expand' index out of range for parameter '%s'",
  //                out_pname.c_str());
  //         }
  //         design_col[r] = cur_design(row_idx, j);
  //       }
  //       const double* d_ptr = design_col.begin();
  //
  //       // Accumulate in-place: out_ptr += coef * design_col, NA/NaN -> 0
  //       for (int r = 0; r < n_trials_local; ++r) {
  //         double v = coef_ptr[r] * d_ptr[r];
  //         if (NumericVector::is_na(v) || std::isnan(v)) {
  //           v = 0.0;
  //         }
  //         out_ptr[r] += v;
  //       }
  //     }
  //   }
  // }

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

