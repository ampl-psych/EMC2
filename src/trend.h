#ifndef dynamic_h
#define dynamic_h

#include "utility_functions.h"
#include "EMC2/userfun.hpp"
#include <Rcpp.h>
#include <unordered_map>
#include <regex>
#include "trend_DBM.h"
using namespace Rcpp;

// Call a user-supplied custom trend kernel
// [[Rcpp::export]]
NumericVector EMC2_call_custom_trend(NumericMatrix trend_pars,
                                     NumericMatrix input,
                                     SEXP funptrSEXP) {
  XPtr<userfun_t> funptr(funptrSEXP);
  if (funptr.get() == nullptr) stop("Null function pointer.");
  userfun_t f = *funptr;
  if (!f) stop("Invalid function pointer.");
  NumericVector res = f(trend_pars, input);
  if (res.size() != input.nrow()) {
    stop("Custom trend function must return a vector of length nrow(input).");
  }
  return res;
}




NumericVector run_delta_rcpp(NumericVector q0, NumericVector alpha, NumericVector covariate) {
  int n = covariate.length();
  NumericVector q(n);
  NumericVector pe(n);
  q[0] = q0[0];
  // forward-looking
  for(int i = 0; i < n-1; i++) {
    if(NumericVector::is_na(covariate[i])) {
      q[i+1] = q[i];
    } else {
      pe[i] = covariate[i]-q[i];
      q[i+1] = q[i] + alpha[i] * pe[i];
    }
  }
  return q;
}

NumericVector run_delta2kernel_rcpp(NumericVector q0, NumericVector alphaFast,
                                    NumericVector propSlow, NumericVector dSwitch,
                                    NumericVector covariate) {
  int n = covariate.length();
  NumericVector q(n);
  NumericVector qFast(n);
  NumericVector qSlow(n);
  NumericVector peFast(n);
  NumericVector peSlow(n);

  q[0] = qFast[0] = qSlow[0] = q0[0];
  NumericVector alphaSlow = propSlow * alphaFast;

  for(int i = 0; i < n-1; i++) {
    if(NumericVector::is_na(covariate[i])) {
      qFast[i+1] = qFast[i];
      qSlow[i+1] = qSlow[i];
      q[i+1] = q[i];
    } else {
      peFast[i] = covariate[i] - qFast[i];
      peSlow[i] = covariate[i] - qSlow[i];
      qFast[i+1] = qFast[i] + alphaFast[i] * peFast[i];
      qSlow[i+1] = qSlow[i] + alphaSlow[i] * peSlow[i];
      if(std::abs(qFast[i+1] - qSlow[i+1]) > dSwitch[i+1]) {
        q[i+1] = qFast[i+1];
      } else {
        q[i+1] = qSlow[i+1];
      }
    }
  }

  return q;
}

NumericVector run_delta2lr_rcpp(NumericVector q0, NumericVector alphaPos,
                                NumericVector alphaNeg,
                                NumericVector covariate) {
  int n = covariate.length();
  NumericVector q(n);
  NumericVector pe(n);
  double alpha;
  q[0] = q0[0];
  for(int i = 0; i < n-1; i++) {
    if(NumericVector::is_na(covariate[i])) {
      q[i+1] = q[i];
    } else {
      pe[i] = covariate[i] - q[i];
      alpha = (pe[i] > 0) ? alphaPos[i] : alphaNeg[i];
      q[i+1] = q[i] + alpha * pe[i];
    }
  }

  return q;
}

// Build expand indices for compress/expand logic from a first-level mask
inline IntegerVector build_expand_idx_rcpp(const LogicalVector& first_level) {
  const int n = first_level.size();
  IntegerVector expand_idx(n);
  int count = 0;
  for (int i = 0; i < n; ++i) {
    if (first_level[i]) ++count;
    expand_idx[i] = count;
  }
  for (int i = 0; i < n; ++i) {
    if (expand_idx[i] == 0) {
      stop("Found rows before first 'at' level within subject. Cannot anchor expansion.");
    }
  }
  return expand_idx;
}


NumericMatrix run_kernel_rcpp(NumericMatrix kernel_pars,
                              String kernel,
                              NumericMatrix input,
                              SEXP funptrSEXP = R_NilValue,
                              LogicalVector first_level_mask = LogicalVector(),
                              bool ffill_na = false) {
  // Kernels accept any number of input columns; apply per column and sum contributions.
  const int n = input.nrow();
  const int p = input.ncol();
  NumericMatrix out;
  if (kernel == "custom") {
    out = NumericMatrix(n, 1);
  } else {
    out = NumericMatrix(n, p);
  }
  out.fill(0.0);

  const bool use_at = (first_level_mask.size() == n);
  IntegerVector expand_idx;
  LogicalVector comp_rows;
  if (use_at) {
    comp_rows = first_level_mask;
    expand_idx = build_expand_idx_rcpp(first_level_mask);
  } else {
    comp_rows = LogicalVector(n, true);
    expand_idx = seq(1, n);
  }

  // Compressed inputs/params if at is used
  NumericMatrix input_comp = submat_rcpp(input, comp_rows);
  NumericMatrix kp_comp = submat_rcpp(kernel_pars, comp_rows);
  const int n_comp = input_comp.nrow();

  // Custom kernel path: take all inputs at once (matrix), exclude rows with any NA, no map
  if (kernel == "custom") {
    if (Rf_isNull(funptrSEXP)) stop("Missing function pointer for custom kernel.");
    // // Row-wise NA check
    // LogicalVector good(n_comp, true);
    // for (int i = 0; i < n_comp; ++i) {
    //   for (int c = 0; c < p; ++c) {
    //     if (NumericVector::is_na(input_comp(i, c))) { good[i] = false; break; }
    //   }
    // }
    // const int n_good = sum(good);
    // NumericVector comp_out(n_comp); // zeros by default
    // if (n_good > 0) {
    //   NumericMatrix in_good(n_good, p);
    //   int rg = 0;
    //   for (int i = 0; i < n_comp; ++i) if (good[i]) {
    //     for (int c = 0; c < p; ++c) in_good(rg, c) = input_comp(i, c);
    //     ++rg;
    //   }
    //   NumericMatrix kp_good = submat_rcpp(kp_comp, good);
    //
    //   NumericVector contrib;
    //   contrib = EMC2_call_custom_trend(kp_good, in_good, funptrSEXP);
    //   if (contrib.size() != n_good) stop("Custom kernel returned wrong length (expected n_good).");
    //   rg = 0;
    //   // for (int i = 0; i < n_comp; ++i) if (good[i]) comp_out[i] = NumericVector::is_na(contrib[rg]) ? 0.0 : contrib[rg++];
    //   // Fill
    //   double last_val = 0.0;
    //   for (int i = 0; i < n_comp; ++i) {
    //     if (good[i]) {
    //       // Take next contrib, replace NA with 0
    //       double val = NumericVector::is_na(contrib[rg]) ? 0.0 : contrib[rg];
    //       comp_out[i] = val;
    //       if (ffill_na) last_val = val;  // remember for forward fill
    //       rg++;
    //     } else if (ffill_na) {
    //       // Forward-fill from last known value
    //       comp_out[i] = last_val;
    //       // SM: This needs a rethink.
    //       // Probably best to just let the user handle NA-values in the custom kernel?
    //       // We cannot know whether there's some sort of delta-rule logic where a backward fill would be needed
    //       // or forward fill logic...
    //     }
    //   }
    // }

    // No NA-check, handle NAs in custom kernel
    NumericVector contrib;
    contrib = EMC2_call_custom_trend(kp_comp, input_comp, funptrSEXP);
    if (contrib.size() != n_comp) stop("Custom kernel returned wrong length (expected n_good).");

    // expand back and return
    for (int i = 0; i < n; ++i) {
      int idx = expand_idx[i] - 1;
      out(i,0) += contrib[idx];
      // out(i,0) += comp_out[idx];
    }
    return out;
  }

  // for convenience, define some regular expressions for matching input kernel
  // against specific classes of kernels
  static const std::string kern_c = kernel.get_cstring();
  static const std::regex delta_kerns_re("^delta(2lr|2kernel)?$");
  static const std::regex beta_binom_kerns_re("^beta_binom(_map)?(_surprise)?$");
  static const std::regex dbm_kerns_re("^dbm(_map)?(_surprise)?$");
  static const std::regex tpm_kerns_re("^tpm(_surprise)?$");
  static const std::regex learning_kerns_re(
      "^(delta(2lr|2kernel)?|"
      "beta_binom(_map)?(_surprise)?|"
      "dbm(_map)?(_surprise)?|"
      "tpm(_surprise)?)$"
  );

  for (int c = 0; c < p; ++c) {
    NumericVector cov_comp = input_comp(_, c);
    NumericVector comp_out(n_comp); // zeros by default

    if(std::regex_match(kern_c, learning_kerns_re)) {
      // Sequential kernels: No pre-filtering of NA, handle NA values in kernel
      if (kernel == "delta") {
        comp_out = run_delta_rcpp(
          // q0, alpha
          kp_comp(_,0), kp_comp(_,1),
          cov_comp
        );
      } else if (kernel == "delta2kernel") {
        comp_out = run_delta2kernel_rcpp(
          // q0, alphaSlow, propSlow, dSwitch
          kp_comp(_,0), kp_comp(_,1), kp_comp(_,2), kp_comp(_,3),
          cov_comp
        );
      } else if (kernel == "delta2lr") {
        comp_out = run_delta2lr_rcpp(
          // q0, alphaPos, alphaNeg
          kp_comp(_,0), kp_comp(_,1), kp_comp(_,2),
          cov_comp);
      } else if (std::regex_match(kern_c, beta_binom_kerns_re)) {
        // beta binomial learning model variants
        bool return_map = false, return_surprise = false;
        if (std::regex_match(kern_c, std::regex("_surprise$"))) {
          return_surprise = true;
        }
        if (std::regex_match(kern_c, std::regex("^beta_binom_map"))) {
          return_map = true;
        }
        comp_out = run_beta_binomial(
          cov_comp,
          // a0, b0, decay, window
          kp_comp(_,0), kp_comp(_,1), kp_comp(_,2), kp_comp(_,3),
          return_map, return_surprise
        );
      } else if (std::regex_match(kern_c, dbm_kerns_re)) {
        // Dynamic Belief Model variants
        bool return_map = false, return_surprise = false;
        if (std::regex_match(kern_c, std::regex("_surprise$"))) {
          return_surprise = true;
        }
        if (std::regex_match(kern_c, std::regex("^dbm_map"))) {
          return_map = true;
        }
        comp_out = run_dbm(
          cov_comp,
          // cp, mu0, s0
          kp_comp(_,0), kp_comp(_,1), kp_comp(_,2),
          return_map, return_surprise
        );
      } else if (std::regex_match(kern_c, tpm_kerns_re)) {
        // Transition Probability Model variants
        const bool return_surprise = (kernel == "tpm_surprise") ? true : false;
        comp_out = run_tpm(
          cov_comp,
          // cp, a0, b0
          kp_comp(_,0), kp_comp(_,1), kp_comp(_,2),
          return_surprise
        );
      }

      if(!ffill_na) {
        // If, for whatever reason, the user wants NA-values to be set to 0, we can still do this.
        LogicalVector good = !is_na(cov_comp);
        // otherwise, apply forward fill
        for (int i = 1; i < n_comp; ++i) {
          if(!good[i]) comp_out[i] = 0;
        }
      }

    } else {

      // non-sequential kernels: Apply NA filtering and ffill logic if needed
      LogicalVector good = !is_na(cov_comp);
      const int n_good = sum(good);
      if (n_good > 0) {
        // Build param subset for good rows (can be 0 columns)
        NumericMatrix kp_good = submat_rcpp(kp_comp, good);
        // Map from comp row index -> position within good subset
        std::vector<int> good_pos(n_comp, -1);
        for (int i = 0, pos = 0; i < n_comp; ++i) if (good[i]) good_pos[i] = pos++;

        // Built-in kernels on good rows only
        if (kernel == "lin_decr") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) comp_out[i] = -cov_comp[i];
        } else if (kernel == "lin_incr") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) comp_out[i] = cov_comp[i];
        } else if (kernel == "exp_decr") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) {
            int r = good_pos[i];
            comp_out[i] = std::exp(-kp_good(r, 0) * cov_comp[i]);
          }
        } else if (kernel == "exp_incr") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) {
            int r = good_pos[i];
            comp_out[i] = 1 - std::exp(-kp_good(r, 0) * cov_comp[i]);
          }
        } else if (kernel == "pow_decr") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) {
            int r = good_pos[i];
            comp_out[i] = std::pow(1 + cov_comp[i], -kp_good(r, 0));
          }
        } else if (kernel == "pow_incr") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) {
            int r = good_pos[i];
            comp_out[i] = 1 - std::pow(1 + cov_comp[i], -kp_good(r, 0));
          }
        } else if (kernel == "poly2") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) {
            int r = good_pos[i];
            comp_out[i] = kp_good(r, 0) * cov_comp[i] + kp_good(r, 1) * std::pow(cov_comp[i], 2);
          }
        } else if (kernel == "poly3") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) {
            int r = good_pos[i];
            comp_out[i] = kp_good(r, 0) * cov_comp[i] + kp_good(r, 1) * std::pow(cov_comp[i], 2) + kp_good(r, 2) * std::pow(cov_comp[i], 3);
          }
        } else if (kernel == "poly4") {
          for (int i = 0; i < n_comp; ++i) if (good[i]) {
            int r = good_pos[i];
            comp_out[i] = kp_good(r, 0) * cov_comp[i] + kp_good(r, 1) * std::pow(cov_comp[i], 2) + kp_good(r, 2) * std::pow(cov_comp[i], 3) + kp_good(r, 3) * std::pow(cov_comp[i], 4);
          }
        } else {
          stop("Unknown kernel type");
        }

        // Fill
        if(ffill_na) {
          // otherwise, apply forward fill
          for (int i = 1; i < n_comp; ++i) {
            if(!good[i]) comp_out[i] = comp_out[i-1];
          }
        }
      }
    }

    // expand comp_out back to full n rows and add
    for (int i = 0; i < n; ++i) {
      int idx = expand_idx[i] - 1; // 0-based
      out(i,c) += comp_out[idx];
    }
  }

  return out;
}

// Now accepts the full parameter matrix `pars_full` so we can use par_input columns as inputs too.
// Passes all inputs (covariates + par_input) to kernel in one call; kernel sums across columns.
// [[Rcpp::export]]
NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param, NumericMatrix trend_pars, NumericMatrix pars_full, bool return_kernel = false) {
  String kernel = as<String>(trend["kernel"]);
  String base = as<String>(trend["base"]);
  // Extract optional custom pointer attribute if present
  SEXP custom_ptr = R_NilValue;
  if (kernel == "custom") {
    custom_ptr = trend.attr("custom_ptr");
  }
  CharacterVector covnames;
  if (trend.containsElementNamed("covariate") && !Rf_isNull(trend["covariate"])) {
    covnames = trend["covariate"];
  } else {
    covnames = CharacterVector(0);
  }
  CharacterVector par_input;
  if (trend.containsElementNamed("par_input") && !Rf_isNull(trend["par_input"])) {
    par_input = trend["par_input"];
  } else {
    par_input = CharacterVector(0);
  }

  bool ffill_na;
  if (trend.containsElementNamed("ffill_na") && !Rf_isNull(trend["ffill_na"])) {
    ffill_na = as<bool>(trend["ffill_na"]);
  } else {
    ffill_na = false;
  }

  // Initialize output vector with zeros
  int n_trials = param.length();
  NumericVector out(n_trials, 0.0);

  int n_base_pars = 0;
  if(base == "lin" || base == "exp_lin" || base == "centered") {
    n_base_pars = 1;
  }
  // Check for covariate maps
  CharacterVector map_names;
  List map_list;
  int n_maps = 0;
  bool has_maps = trend.containsElementNamed("map") && !Rf_isNull(trend["map"]);
  if(has_maps) {
    map_list = trend["map"];
    map_names = map_list.names();
    n_maps = map_names.size();
    n_base_pars = n_maps*n_base_pars; // find n maps
  }

  // extract kernel parameters
  int n_trend_pars = trend_pars.ncol();
  NumericMatrix kernel_pars(n_trials, n_trend_pars- n_base_pars);
  for (int j = n_base_pars; j < n_trend_pars; j++) {
    kernel_pars(_, j - n_base_pars) = trend_pars(_, j);
  }

  // Build input matrix from covariates and par_input columns
  int n_cov = covnames.size();
  // Keep only par_input columns that actually exist in pars_full (if provided)
  CharacterVector pars_full_names = colnames(pars_full);
  std::vector<std::string> par_in_keep;
  for (int i = 0; i < par_input.size(); i++) {
    std::string nm = Rcpp::as<std::string>(par_input[i]);
    bool found = false;
    for (int j = 0; j < pars_full_names.size(); j++) {
      if (nm == Rcpp::as<std::string>(pars_full_names[j])) { found = true; break; }
    }
    if (found) par_in_keep.push_back(nm);
  }
  int n_par_in = (int)par_in_keep.size();

  int n_inputs = n_cov + n_par_in;
  // If no inputs provided, keep out as zeros (base will handle accordingly)
  if (n_inputs > 0) {
    NumericMatrix input_all(n_trials, n_inputs);
    // Fill covariate columns first
    for (int i = 0; i < n_cov; i++) {
      String cur_cov = covnames[i];
      NumericVector covariate = as<NumericVector>(data[cur_cov]);
      input_all(_, i) = covariate;
    }
    // Then par_input columns
    if (n_par_in > 0) {
      CharacterVector par_in_keep_cv(n_par_in);
      for (int i = 0; i < n_par_in; i++) par_in_keep_cv[i] = par_in_keep[i];
      NumericMatrix pin = submat_rcpp_col_by_names(pars_full, par_in_keep_cv);
      for (int j = 0; j < n_par_in; j++) {
        input_all(_, n_cov + j) = pin(_, j);
      }
    }

    // Build optional first-level mask and call kernel (handles compression/expansion internally)
    LogicalVector first_level;
    if (trend.containsElementNamed("at") && !Rf_isNull(trend["at"])) {
      String at_name = trend["at"];
      SEXP at_col = data[at_name];
      if (!Rf_inherits(at_col, "factor")) stop("'at' column must be a factor");
      IntegerVector f = as<IntegerVector>(at_col);
      first_level = (f == 1);
    } else {
      first_level = LogicalVector(0);
    }

    // Run kernel
    NumericMatrix kernel_out = run_kernel_rcpp(kernel_pars, kernel, input_all, custom_ptr, first_level, ffill_na);
    if (return_kernel) return kernel_out;

    int n_rows = kernel_out.nrow();
    int n_cols = kernel_out.ncol();

    // Outer loop over maps
    List covariate_maps = data.attr("covariate_maps");
    int n_loops = has_maps ? n_maps : 1;
    for (int map_n = 0; map_n < n_loops; ++map_n) {

      NumericMatrix covariate_map;
      if (has_maps) {
        std::string map_name = as<std::string>(map_names[map_n]);
        if (Rf_isNull(covariate_maps[map_name])) {
          stop("No map named '%s' found in covariate_maps", map_name.c_str());
        }
        covariate_map = as<NumericMatrix>(covariate_maps[map_name]);
      }

      // Loop over columns
      for (int c = 0; c < n_cols; ++c) {
        // Loop over rows
        for (int r = 0; r < n_rows; ++r) {
          double contrib = kernel_out(r, c);

          // Multiply by covariate map if it exists
          if (has_maps) contrib *= covariate_map(r, c);

          // Apply base parameter
          if (base == "lin" || base == "exp_lin") {
            contrib *= trend_pars(r, map_n);
          } else if (base == "centered") {
            contrib = (contrib - 0.5) * trend_pars(r, map_n);
          } // "add" leaves contrib unchanged

          // Accumulate into output
          out(r) += contrib;
        }
      }
    }
  }

  // Add parameter to obtain final summed input
  if(base == "lin") {
    out += param;
  }
  else if(base == "exp_lin") {
    out += exp(param);
  }
  else if(base == "centered") {
    out += param;
  }
  else if(base == "add") {
    out += param;
  }
  return out;
}

// A few unneccessary loops in here, but seems reasonably efficient
inline NumericMatrix prep_trend_phase(DataFrame data, List trend, NumericMatrix pars, String phase) {
  CharacterVector trend_names = trend.names();
  CharacterVector par_names = colnames(pars);
  CharacterVector all_remove;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    std::string cur_ph = Rcpp::as<std::string>(cur_trend["phase"]);
    std::string phase_s = Rcpp::as<std::string>(Rcpp::wrap(phase));
    if (cur_ph != phase_s) continue;
    String par = trend_names[i];
    CharacterVector trend_pnames = cur_trend["trend_pnames"];
    all_remove = c_add_charvectors(all_remove, trend_pnames);
    LogicalVector par_idx = contains(par_names, par);
    NumericVector param = as<NumericVector>(submat_rcpp_col(pars, par_idx));
    NumericMatrix trend_pars = submat_rcpp_col_by_names(pars, trend_pnames);
    pars(_, as<int>(which_rcpp(par_idx))) = run_trend_rcpp(data, cur_trend, param, trend_pars, pars);
  }
  all_remove = unique(all_remove);
  if (all_remove.size() > 0) {
    CharacterVector pnames = colnames(pars);
    LogicalVector keep = !contains_multiple(pnames, all_remove);
    pars = submat_rcpp_col(pars, keep);
  }
  return(pars);
}

inline NumericMatrix prep_trend_phase_with_pars(DataFrame data, List trend, NumericMatrix pars, String phase, const NumericMatrix& trend_pars_override) {
  CharacterVector trend_names = trend.names();
  CharacterVector par_names = colnames(pars);
  CharacterVector all_remove;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    // Simple cast to std::string for comparison
    std::string cur_ph = Rcpp::as<std::string>(cur_trend["phase"]);
    std::string phase_s = Rcpp::as<std::string>(Rcpp::wrap(phase));
    if (cur_ph != phase_s) continue;
    String par = trend_names[i];
    CharacterVector trend_pnames = cur_trend["trend_pnames"];
    all_remove = c_add_charvectors(all_remove, trend_pnames);
    LogicalVector par_idx = contains(par_names, par);
    NumericVector param = as<NumericVector>(submat_rcpp_col(pars, par_idx));
    NumericMatrix cur_tp = submat_rcpp_col_by_names(trend_pars_override, trend_pnames);
    pars(_, as<int>(which_rcpp(par_idx))) = run_trend_rcpp(data, cur_trend, param, cur_tp, pars);
  }
  all_remove = unique(all_remove);
  if (all_remove.size() > 0) {
    CharacterVector pnames = colnames(pars);
    LogicalVector keep = !contains_multiple(pnames, all_remove);
    pars = submat_rcpp_col(pars, keep);
  }
  return(pars);
}

// ---- Trend helpers for mapping pipeline ----

// Collect all unique trend parameter names across trend list entries
inline CharacterVector collect_trend_param_names(const List& trend) {
  CharacterVector trend_pnames;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    CharacterVector cur_names = cur_trend["trend_pnames"];
    trend_pnames = c_add_charvectors(trend_pnames, cur_names);
  }
  return unique(trend_pnames);
}

inline CharacterVector collect_trend_param_names_phase(const List& trend, const std::string& phase) {
  CharacterVector trend_pnames;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    std::string ph = Rcpp::as<std::string>(cur_trend["phase"]);
    if (ph != phase) continue;
    CharacterVector cur_names = cur_trend["trend_pnames"];
    trend_pnames = c_add_charvectors(trend_pnames, cur_names);
  }
  return unique(trend_pnames);
}

// Build per-trial columns for trend parameters from designs and p_vector, and apply transforms
inline NumericMatrix build_trend_columns_from_design(NumericVector p_vector,
                                                     CharacterVector p_types,
                                                     List designs,
                                                     int n_trials,
                                                     const List& trend,
                                                     const std::vector<TransformSpec>& full_specs) {
  CharacterVector trend_pnames = collect_trend_param_names(trend);
  if (trend_pnames.size() == 0) {
    return NumericMatrix(n_trials, 0); // empty
  }

  // Map parameter name -> index in p_types/designs
  std::unordered_map<std::string,int> name_to_idx;
  for (int i = 0; i < p_types.size(); ++i) {
    name_to_idx[ Rcpp::as<std::string>(p_types[i]) ] = i;
  }

  // Allocate output matrix with columns strictly in trend_pnames order
  NumericMatrix trend_pars(n_trials, trend_pnames.size());
  colnames(trend_pars) = trend_pnames;

  for (int c = 0; c < trend_pnames.size(); ++c) {
    std::string pname = Rcpp::as<std::string>(trend_pnames[c]);
    auto it = name_to_idx.find(pname);
    int idx = it->second;
    NumericMatrix cur_design = designs[idx];
    CharacterVector cur_names = colnames(cur_design);

    // Accumulate p_vector * design columns
    NumericVector acc(n_trials, 0.0);
    for (int j = 0; j < cur_design.ncol(); ++j) {
      String cur_name(cur_names[j]);
      NumericVector tmp = p_vector[cur_name] * cur_design(_, j);
      LogicalVector bad = is_na(tmp) | is_nan(tmp);
      tmp[bad] = 0;
      acc = acc + tmp;
    }
    trend_pars(_, c) = acc;
  }

  // Transform trend parameter columns using precomputed specs
  std::vector<TransformSpec> t_specs = make_transform_specs_from_full(trend_pars, p_types, full_specs);
  trend_pars = c_do_transform(trend_pars, t_specs);
  return trend_pars;
}

// Apply premap trend for a single parameter vector if a trend is defined
inline NumericVector apply_premap_trends(const DataFrame& data,
                                         const List& trend,
                                         const CharacterVector& trend_names,
                                         const String& param_name,
                                         NumericVector param_values,
                                         const NumericMatrix& trend_pars,
                                         const NumericVector& p_vector) {
  // Apply all premap trend entries for this param_name sequentially, in order
  NumericVector result = clone(param_values);
  for (int i = 0; i < trend.size(); ++i) {
    if (trend_names[i] == param_name) {
      List cur_trend = trend[i];
      std::string ph = Rcpp::as<std::string>(cur_trend["phase"]);
      if (ph != "premap") continue;
      CharacterVector cur_trend_pnames = cur_trend["trend_pnames"];
      NumericMatrix cur_trend_pars = submat_rcpp_col_by_names(trend_pars, cur_trend_pnames);
      // Build a small pars_full matrix for par_input: replicate scalars from p_vector across trials
      CharacterVector par_input;
      if (cur_trend.containsElementNamed("par_input") && !Rf_isNull(cur_trend["par_input"])) {
        par_input = cur_trend["par_input"];
      } else {
        par_input = CharacterVector(0);
      }
      NumericMatrix pars_full(result.size(), par_input.size());
      if (par_input.size() > 0) {
        colnames(pars_full) = par_input;
        CharacterVector pnames = p_vector.names();
        for (int c = 0; c < par_input.size(); ++c) {
          std::string nm = Rcpp::as<std::string>(par_input[c]);
          // find scalar in p_vector by name
          int idx = -1;
          for (int k = 0; k < pnames.size(); ++k) {
            if (nm == Rcpp::as<std::string>(pnames[k])) { idx = k; break; }
          }
          double val = (idx >= 0) ? (double)p_vector[idx] : NA_REAL;
          for (int r = 0; r < pars_full.nrow(); ++r) pars_full(r, c) = val;
        }
      }
      result = run_trend_rcpp(data, cur_trend, result, cur_trend_pars, pars_full);
    }
  }
  return result;
}

// Fill columns for trend parameters into the output matrix by name (for pretransform case)
inline void fill_trend_columns_for_pretransform(NumericMatrix& pars,
                                                const CharacterVector& p_types,
                                                const NumericMatrix& trend_pars) {
  if (trend_pars.ncol() == 0) return;
  CharacterVector tnames = colnames(trend_pars);
  // Build a name->col map for trend_pars
  std::unordered_map<std::string,int> tmap;
  for (int c = 0; c < tnames.size(); ++c) {
    tmap[Rcpp::as<std::string>(tnames[c])] = c;
  }
  for (int i = 0; i < p_types.size(); ++i) {
    std::string pname = Rcpp::as<std::string>(p_types[i]);
    auto it = tmap.find(pname);
    if (it != tmap.end()) {
      pars(_, i) = trend_pars(_, it->second);
    }
  }
}

#endif


// Rprintf("Subject %d", subject_number);
// Rprintf("Data shape (nrow) %d", n_trials);
// int n_show = std::min(5, map_in.nrow());  // show first 5 rows
// Rcout << "First " << n_show << " rows of map_in:\n";
//
// for (int i = 0; i < n_show; i++) {
//   for (int j = 0; j < map_in.ncol(); j++) {
//     Rcout << map_in(i, j) << "\t";
//   }
//   Rcout << "\n";
// }
// Rcout << std::endl;

// int subject_number = as<IntegerVector>(data["subjects"])[0];
// NumericMatrix map_in(n_trials, n_cov);
// int data_r = 0;
// for(int r = 0; r < map_provided.nrow(); r++) {
//   if(map_provided(r,0) == subject_number) {
//     for(int c = 0; c < n_cov; c++) {
//       map_in(data_r, c) = map_provided(r, c+1);
//     }
//     data_r++;
//   }
// }

// NumericMatrix covariate_map;
// if(has_map) {
//   // find all maps to be applied
//   List map_list = trend["map"];
//   CharacterVector map_names = map_list.names();
//   int n_maps = map_names.size();
//   if (n_maps == 0) {
//     stop("trend$map has no named elements");
//   }
//
//   // Loop over covariate_maps to be applied to this kernel
//   List covariate_maps = data.attr("covariate_maps");
//   for(int map_n = 0; map_n < n_maps; map_n++) {
//     std::string map_name = as<std::string>(map_names[map_n]);
//     if (Rf_isNull(covariate_maps[map_name])) {
//       stop("No map named '%s' found in covariate_maps", map_name.c_str());
//     }
//     covariate_map = as<NumericMatrix>(covariate_maps[map_name]);
//     for (int c = 0; c < kernel_out.ncol(); ++c) {
//       if(base == "lin" || base == "exp_lin") {
//         out += kernel_out(_, c) * covariate_map(_, c) * trend_pars(_,c);
//       } else if(base == "centered") {
//         out += ((kernel_out(_, c) * covariate_map(_, c)) - 0.5) * trend_pars(_,c);
//       } else if(base == "add") {
//         out += kernel_out(_, c) * covariate_map(_, c);
//       }
//     }
//   }
// } else {
//   // No maps provided
//   for (int c = 0; c < kernel_out.ncol(); ++c) {
//     if(base == "lin" || base == "exp_lin") {
//       out += kernel_out(_, c) * trend_pars(_,c);
//     } else if(base == "centered") {
//       out += (kernel_out(_, c) - 0.5) * trend_pars(_,c);
//     } else if(base == "add") {
//       out += kernel_out(_, c);
//     }
//   }
// }
