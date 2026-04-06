#pragma once

#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
// #include <RcppArmadillo.h>
using Rcpp::_;
// using arma::mat;
// using arma::vec;

struct DesignEntry {
  bool valid;                    // does this design have a usable mapping?
  bool skip_self_intercept;      // from design_is_self_intercept
  int out_idx;                   // column index in base to write into
  std::vector<int> coef_idx;     // length K; -1 for unused
  bool uses_self;                // does any coef refer to out_idx?
};

// View-based ParamTable: one base matrix + active column indices
// Invariants:
// - base contains all parameters (never shrinks)
// - active_cols is a view into base (ordering may change)
// - name_to_base_idx maps ALL base columns, active or not
struct ParamTable {
  Rcpp::NumericMatrix base;     // underlying matrix (all columns)

  std::vector<DesignEntry> design_plan;   // cached once

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

  void init_design_plan(const Rcpp::List& designs) {
    using namespace Rcpp;
    using std::string;

    const int n_params = designs.size();
    design_plan.clear();
    design_plan.resize(n_params);

    CharacterVector design_names = designs.names();
    if (design_names.size() != n_params) {
      stop("ParamTable::init_design_plan: designs must be a named list");
    }

    const int T = base.nrow();

    for (int i = 0; i < n_params; ++i) {
      DesignEntry& entry = design_plan[i];
      entry.valid = false;

      if (designs[i] == R_NilValue) continue;

      // self-intercept-only designs can be marked as skippable
      bool skip_self = !design_is_self_intercept.empty() &&
        i < (int)design_is_self_intercept.size() &&
        design_is_self_intercept[i];
      entry.skip_self_intercept = skip_self;

      const string out_name = as<string>(design_names[i]);
      auto out_it = name_to_base_idx.find(out_name);
      if (out_it == name_to_base_idx.end()) {
        // no matching output column; nothing to do
        continue;
      }
      entry.out_idx = out_it->second;

      NumericMatrix design = designs[i];
      if (design.nrow() != T) {
        stop("ParamTable::init_design_plan: design for '%s' must have n_trials rows",
             out_name.c_str());
      }

      const int K = design.ncol();
      CharacterVector coef_names = colnames(design);

      entry.coef_idx.assign(K, -1);
      entry.uses_self = false;

      for (int j = 0; j < K; ++j) {
        string coef_name = as<string>(coef_names[j]);
        auto it = name_to_base_idx.find(coef_name);
        if (it == name_to_base_idx.end()) continue;

        int cidx = it->second;
        entry.coef_idx[j] = cidx;
        if (cidx == entry.out_idx) entry.uses_self = true;
      }

      entry.valid = true;
    }
  }

  Rcpp::NumericMatrix materialize_by_param_names(const Rcpp::CharacterVector& param_names) const {
    using namespace Rcpp;
    const int k = param_names.size();
    NumericMatrix out(n_trials, k);
    CharacterVector out_names(k);

    for (int j = 0; j < k; ++j) {
      std::string nm = as<std::string>(param_names[j]);
      int base_idx = base_index_for(nm);  // throws if unknown

      // Ensure this column is active (if you care about active_cols)
      // bool is_active = false;
      // for (int a : active_cols) {
      //   if (a == base_idx) {
      //     is_active = true;
      //     break;
      //   }
      // }
      // if (!is_active) {
      //   stop("ParamTable::materialize_by_param_names: parameter '%s' not active", nm.c_str());
      // }

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

    // If no include_param or wrong length, include all
    LogicalVector use =
      (include_param.size() == n_params)
      ? include_param
    : LogicalVector(n_params, true);

    // Lazy initialisation of the plan
    if (design_plan.size() != (size_t)n_params) {
      init_design_plan(designs);
    }

    const int T = n_trials;

    for (int i = 0; i < n_params; ++i) {
      if (!use[i]) continue;
      if (designs[i] == R_NilValue) continue;

      DesignEntry& entry = design_plan[i];
      if (!entry.valid) continue;

      // skip self-intercept-only designs
      if (entry.skip_self_intercept) {
        continue;
      }

      const int out_idx = entry.out_idx;

      NumericMatrix design = designs[i];   // still a light handle
      const int K = design.ncol();

      // Preserve self column if needed
      std::vector<double> self_copy;
      double* out = &base(0, out_idx);

      if (entry.uses_self) {
        self_copy.resize(T);
        const double* src = out;
        std::copy(src, src + T, self_copy.begin());
      }

      // Clear output
      std::fill(out, out + T, 0.0);

      // Main accumulation loop (unchanged logic)
      for (int j = 0; j < K; ++j) {
        int cidx = entry.coef_idx[j];
        if (cidx < 0) continue;

        const double* coef =
          (entry.uses_self && cidx == out_idx)
          ? self_copy.data()
            : &base(0, cidx);

        const double* d = &design(0, j);

        for (int r = 0; r < T; ++r) {
          double v = coef[r] * d[r];
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

  // Or pass a matrix and the row to fill. Should be faster due to the pre-defined
  // pm_col_to_base_idx vector
  void fill_from_particle_row(const Rcpp::NumericMatrix& particles,
                              int row,
                              const std::vector<int>& pm_col_to_base_idx)
  {
    const int ncols = particles.ncol();
    const int n_trials = this->n_trials; // from ParamTable

    // Zero everything once per particle
    reset_base_to_zero();

    for (int j = 0; j < ncols; ++j) {
      int base_idx = pm_col_to_base_idx[j];
      if (base_idx < 0) continue;   // particle has a param we don't use
      double val = particles(row, j);
      double* col = &base(0, base_idx);
      for (int r = 0; r < n_trials; ++r) {
        col[r] = val;
      }
    }
  }
};



std::unordered_set<std::string> param_names_excluding(const ParamTable& pt,
                                                      std::initializer_list<const std::unordered_set<std::string>*> excludes);

Rcpp::CharacterVector names_excluding(const Rcpp::CharacterVector& names,
                                      std::initializer_list<const std::unordered_set<std::string>*> excludes);

Rcpp::NumericMatrix add_constants_columns(Rcpp::NumericMatrix p_matrix,
                                          Rcpp::NumericVector constants);
