#pragma once
#ifndef PARAM_TABLE_H
#define PARAM_TABLE_H

#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include "Mat.h"
using Rcpp::_;


// ---------------------------------------------------------------------------
// DesignEntry
// ---------------------------------------------------------------------------
struct DesignEntry {
  bool valid                = false;
  bool skip_self_intercept  = false;
  int  out_idx              = -1;
  std::vector<int> coef_idx;
  bool uses_self            = false;
  bool dm_is_constant       = false;  // true if every column of the design matrix is constant

  // Raw design matrix data — pointer into the original Rcpp::NumericMatrix,
  // which is owned externally and must outlive this entry (safe: R keeps it alive).
  const double* dm_ptr  = nullptr;
  int           dm_nrow = 0;
  int           dm_ncol = 0;
  bool          dm_null = true;
};

// ---------------------------------------------------------------------------
// ParamTable: one Mat base + active column indices
// ---------------------------------------------------------------------------
struct ParamTable {
  Mat base;                                        // underlying matrix (all columns)
  std::vector<DesignEntry> design_plan;            // cached once per setup

  std::vector<std::string> base_names;             // colnames for base (plain strings)
  int n_trials = 0;

  std::vector<int>  active_cols;
  std::vector<char> is_active;
  std::vector<bool> col_is_constant;   // true if all rows of this column are identical

  std::unordered_map<std::string, int> name_to_base_idx;

  ParamTable() = default;

  ParamTable(Mat base_, std::vector<std::string> names_)
    : base(std::move(base_)), base_names(std::move(names_))
  {
    n_trials = base.nrow;
    const int p = base.ncol;
    active_cols.resize(p);
    std::iota(active_cols.begin(), active_cols.end(), 0);
    col_is_constant.assign(p, false);  // tracks if all values are constants
    rebuild_name_map();
  }

  // Construct from Rcpp types — call this at the R boundary, before any parallel region
  static ParamTable from_rcpp(const Rcpp::NumericMatrix& rcpp_base,
                              const Rcpp::CharacterVector& rcpp_names)
  {
    Mat m = Mat::from_rcpp(rcpp_base);
    std::vector<std::string> names(rcpp_names.size());
    for (int i = 0; i < (int)rcpp_names.size(); ++i)
      names[i] = Rcpp::as<std::string>(rcpp_names[i]);
    return ParamTable(std::move(m), std::move(names));
  }

  ParamTable deep_copy() const {
    ParamTable pt;
    pt.base                     = base.clone();
    pt.base_names               = base_names;
    pt.n_trials                 = n_trials;
    pt.active_cols              = active_cols;
    pt.is_active                = is_active;
    pt.design_plan              = design_plan;
    pt.col_is_constant          = col_is_constant;
    pt.rebuild_name_map();
    return pt;
  }

  void rebuild_name_map() {
    name_to_base_idx.clear();
    name_to_base_idx.reserve(base_names.size());
    for (int j = 0; j < (int)base_names.size(); ++j)
      name_to_base_idx[base_names[j]] = j;
  }

  int base_index_for(const std::string& nm) const {
    auto it = name_to_base_idx.find(nm);
    if (it == name_to_base_idx.end())
      Rcpp::stop("ParamTable: Unknown parameter '%s'", nm.c_str());
    return it->second;
  }

  int n_params() const { return (int)active_cols.size(); }

  const std::string& name_at(int j) const { return base_names[active_cols[j]]; }

  bool is_active_base_idx(int base_idx) const {
    return std::find(active_cols.begin(), active_cols.end(), base_idx) != active_cols.end();
  }

  const double* column_by_name_ptr(const std::string& name) const {
    return base.colptr(base_index_for(name));
  }

  // Assign into existing column (writes into base)
  void set_column(const std::string& nm, const double* src) {
    int idx = base_index_for(nm);
    double* dst = base.colptr(idx);
    std::copy(src, src + n_trials, dst);
  }

  // For code that still needs a Rcpp view (non-hot paths only)
  Rcpp::NumericVector column_by_name_rcpp(const std::string& nm) const {
    int idx = base_index_for(nm);
    const double* col = base.colptr(idx);
    return Rcpp::NumericVector(col, col + n_trials);
  }

  // ---------------------------------------------------------------------------
  // fill_from_particle_row — raw pointer, OpenMP-safe
  // ---------------------------------------------------------------------------
  void fill_from_particle_row(const double*           particles_ptr,
                              int                     n_particle_cols,
                              int                     row,
                              int                     n_particle_rows,
                              const std::vector<int>& pm_col_to_base_idx)
  {
    reset_base_to_zero();
    for (int j = 0; j < n_particle_cols; ++j) {
      int base_idx = pm_col_to_base_idx[j];
      if (base_idx < 0) continue;
      const double val = particles_ptr[j * n_particle_rows + row];
      double* col = base.colptr(base_idx);
#pragma omp simd
      for (int r = 0; r < n_trials; ++r) col[r] = val;
    }
  }

  // Legacy Rcpp overload — used in serial calc_ll path
  void fill_from_particle_row(const Rcpp::NumericMatrix& particles,
                              int row,
                              const std::vector<int>& pm_col_to_base_idx)
  {
    const double* ptr = particles.begin();
    fill_from_particle_row(ptr, particles.ncol(), row, particles.nrow(), pm_col_to_base_idx);
  }

  // ---------------------------------------------------------------------------
  // fill_from_p_vector
  // ---------------------------------------------------------------------------
  void fill_from_p_vector(const Rcpp::NumericVector& p_vector) {
    Rcpp::CharacterVector pv_names = p_vector.names();
    const int n_names = pv_names.size();
    for (int i = 0; i < n_names; ++i) {
      std::string nm = Rcpp::as<std::string>(pv_names[i]);
      auto it = name_to_base_idx.find(nm);
      if (it == name_to_base_idx.end()) continue;
      double val = p_vector[i];
      double* col = base.colptr(it->second);
#pragma omp simd
      for (int r = 0; r < n_trials; ++r)
        col[r] = val;
    }
  }

  // ---------------------------------------------------------------------------
  // reset_base_to_zero
  // ---------------------------------------------------------------------------
  void reset_base_to_zero() {
    const int N = n_trials * base.ncol;
    double* d = base.data.data();
#pragma omp simd
    for (int i = 0; i < N; ++i)
      d[i] = 0.0;
  }

  // ---------------------------------------------------------------------------
  // map_from_designs — hot path, no Rcpp types touched after init
  // ---------------------------------------------------------------------------
  void map_from_designs(const std::vector<bool>& use)
  {
    const int n = (int)design_plan.size();
    if (n == 0) return;

    // Lazy init of design_plan is not safe inside a parallel region —
    // caller must ensure init_design_plan() was called before entering parallel.
    // (calc_ll_multithreaded calls ctx.param_table.init_design_plan(designs) before the loop)

    const int T = n_trials;
    std::vector<double> self_copy(T);

    for (int i = 0; i < n; ++i) {
      if (!use[i])       continue;

      DesignEntry& entry = design_plan[i];
      if (!entry.valid || entry.skip_self_intercept) continue;
      if (entry.dm_null)    continue;

      const int    out_idx = entry.out_idx;
      const int    K       = entry.dm_ncol;
      double*      out     = base.colptr(out_idx);

      // One of the input names equals the output name -- need to copy the value
      if(entry.uses_self) std::copy(out, out + T, self_copy.begin());

      std::fill(out, out + T, 0.0);

      for (int j = 0; j < K; ++j) {
        int cidx = entry.coef_idx[j];
        if (cidx < 0) continue;

        //
        const double* parameter_value = (entry.uses_self && cidx == out_idx) ? self_copy.data() : base.colptr(cidx);
        const double* design_coef = entry.dm_ptr + j * entry.dm_nrow;

        if (col_is_constant[cidx]) {
          // parameter values are uniform across all rows — scalar fast path
          const double param = parameter_value[0];
          if (param == 0.0) continue;
          if (std::isfinite(param)) {
            // finite values of param guaranteed, so no guard needed
#pragma omp simd
            for (int r = 0; r < T; ++r)
              out[r] += param * design_coef[r];
          } else {
            // ±Inf scalar: guard against 0 * Inf -
            // should be rare, only when parameters are set constant to qnorm(0) or log(0) and
            // the design formula has p~0
#pragma omp simd // might still be used on some compilers?
            for (int r = 0; r < T; ++r) {
              out[r] += (design_coef[r] != 0.0) ? param * design_coef[r] : 0.0;
            }
          }
        } else {
#pragma omp simd
          for (int r = 0; r < T; ++r) {
            // ±Inf parameter value: guard against 0 * Inf -
            // should be rare, only when parameters are set constant to qnorm(0) or log(0) and
            // the design formula has p~0
            out[r] += (design_coef[r] != 0.0) ? parameter_value[r] * design_coef[r] : 0.0;
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // materialize (non-hot, returns Rcpp matrix for R interface)
  // ---------------------------------------------------------------------------
  Rcpp::NumericMatrix materialize() const {
    const int p = (int)active_cols.size();
    Rcpp::NumericMatrix out(n_trials, p);
    Rcpp::CharacterVector out_names(p);
    for (int j = 0; j < p; ++j) {
      int base_j = active_cols[j];
      double* dst = &out(0, j);
      const double* src = base.colptr(base_j);
      std::copy(src, src + n_trials, dst);
      out_names[j] = base_names[base_j];
    }
    Rcpp::colnames(out) = out_names;
    return out;
  }

  Rcpp::NumericMatrix materialize_by_param_names(const Rcpp::CharacterVector& param_names) const {
    const int k = param_names.size();
    Rcpp::NumericMatrix out(n_trials, k);
    Rcpp::CharacterVector out_names(k);
    for (int j = 0; j < k; ++j) {
      std::string nm = Rcpp::as<std::string>(param_names[j]);
      int base_idx = base_index_for(nm);
      double* dst = &out(0, j);
      const double* src = base.colptr(base_idx);
      std::copy(src, src + n_trials, dst);
      out_names[j] = base_names[base_idx];
    }
    Rcpp::colnames(out) = out_names;
    return out;
  }

  // ---------------------------------------------------------------------------
  // Static constructors
  // ---------------------------------------------------------------------------
  static ParamTable from_p_types(int n_trials,
                                 const Rcpp::CharacterVector& p_types)
  {
    const int p = p_types.size();
    std::vector<std::string> names(p);
    for (int i = 0; i < p; ++i)
      names[i] = Rcpp::as<std::string>(p_types[i]);
    return ParamTable(Mat(n_trials, p), std::move(names));
  }

  static ParamTable from_p_vector_and_designs(const Rcpp::NumericVector& p_vector,
                                              const Rcpp::List& designs,
                                              int n_trials)
  {
    using namespace Rcpp;
    using std::string;

    // --- Build name list ---
    std::vector<string> names_vec;
    std::unordered_set<string> seen;
    names_vec.reserve(p_vector.size() + designs.size() * 2);

    CharacterVector pv_names     = p_vector.names();
    CharacterVector design_names = designs.names();

    for (int i = 0; i < pv_names.size(); ++i) {
      string nm = as<string>(pv_names[i]);
      if (seen.insert(nm).second) names_vec.push_back(nm);
    }
    for (int i = 0; i < design_names.size(); ++i) {
      string nm = as<string>(design_names[i]);
      if (seen.insert(nm).second) names_vec.push_back(nm);
    }
    for (int i = 0; i < designs.size(); ++i) {
      if (designs[i] == R_NilValue) continue;
      NumericMatrix dm = designs[i];
      CharacterVector cn = colnames(dm);
      for (int j = 0; j < cn.size(); ++j) {
        string nm = as<string>(cn[j]);
        if (seen.insert(nm).second) names_vec.push_back(nm);
      }
    }

    // --- Build base matrix ---
    const int p = (int)names_vec.size();
    Mat base(n_trials, p);

    std::unordered_map<string, double> pval;
    pval.reserve(p_vector.size());
    for (int i = 0; i < p_vector.size(); ++i)
      pval[as<string>(pv_names[i])] = p_vector[i];

    for (int j = 0; j < p; ++j) {
      auto it = pval.find(names_vec[j]);
      if (it != pval.end()) {
        double val = it->second;
        double* col = base.colptr(j);
        for (int r = 0; r < n_trials; ++r) col[r] = val;
      }
    }

    // --- Construct ParamTable (builds name_to_base_idx) ---
    ParamTable pt(std::move(base), names_vec);
    pt.n_trials = n_trials;
    // --- Build design_plan ---
    const int n_designs = designs.size();
    pt.design_plan.clear();
    pt.design_plan.resize(n_designs);

    for (int i = 0; i < n_designs; ++i) {
      DesignEntry& entry = pt.design_plan[i];
      entry.valid = false;

      if (designs[i] == R_NilValue) continue;

      NumericMatrix   dm       = designs[i];
      const string    out_name = as<string>(design_names[i]);

      auto out_it = pt.name_to_base_idx.find(out_name);
      if (out_it == pt.name_to_base_idx.end()) continue;
      entry.out_idx = out_it->second;

      if (dm.nrow() != n_trials)
        stop("ParamTable: design for '%s' must have n_trials rows", out_name.c_str());

      // --- Raw pointer fields ---
      entry.dm_ptr  = dm.begin();
      entry.dm_nrow = dm.nrow();
      entry.dm_ncol = dm.ncol();
      entry.dm_null = false;

      // --- dm_is_constant: every column has identical rows ---
      entry.dm_is_constant = true;
      for (int j = 0; j < entry.dm_ncol && entry.dm_is_constant; ++j) {
        const double* col = entry.dm_ptr + j * entry.dm_nrow;
        const double  v0  = col[0];
        for (int r = 1; r < entry.dm_nrow; ++r)
          if (col[r] != v0) { entry.dm_is_constant = false; break; }
      }

      // --- skip_self_intercept: single all-ones self-column — mapping is a no-op ---
      entry.skip_self_intercept = false;
      if (entry.dm_ncol == 1) {
        CharacterVector cn = colnames(dm);
        if (as<string>(cn[0]) == out_name) {
          const double* col = entry.dm_ptr;
          bool all_ones = true;
          for (int r = 0; r < entry.dm_nrow; ++r)
            if (col[r] != 1.0) { all_ones = false; break; }
          entry.skip_self_intercept = all_ones;
        }
      }

      // --- Coefficient index mapping ---
      const int K = entry.dm_ncol;
      CharacterVector coef_names = colnames(dm);
      entry.coef_idx.assign(K, -1);
      entry.uses_self = false;

      for (int j = 0; j < K; ++j) {
        string coef_name = as<string>(coef_names[j]);
        auto it = pt.name_to_base_idx.find(coef_name);
        if (it == pt.name_to_base_idx.end()) continue;
        int cidx = it->second;
        entry.coef_idx[j] = cidx;
        if (cidx == entry.out_idx) entry.uses_self = true;
      }

      entry.valid = true;
    }

    return pt;
  }
};

// ---------------------------------------------------------------------------
// Free helpers
// ---------------------------------------------------------------------------
std::unordered_set<std::string> param_names_excluding(
    const ParamTable& pt,
    std::initializer_list<const std::unordered_set<std::string>*> excludes);

Rcpp::CharacterVector names_excluding(
    const Rcpp::CharacterVector& names,
    std::initializer_list<const std::unordered_set<std::string>*> excludes);

Rcpp::NumericMatrix add_constants_columns(Rcpp::NumericMatrix p_matrix,
                                          Rcpp::NumericVector constants);

#endif
