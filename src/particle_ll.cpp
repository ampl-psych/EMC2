#include <Rcpp.h>
#include <unordered_map>

// Utilities first — no dependencies on model types
#include "utility_functions.h"
#include "transform_utils.h"
#include "ParamTable.h"
#include "trend.h"
#include "TrendEngine.h"
#include "math_utils.h"

// Model headers — each includes RaceSpec.h themselves
#include "model_LBA.h"
#include "model_lnr.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"

// RaceSetup last — references functions defined in model headers above
#include "RaceSetup.h"
using namespace Rcpp;


// cache for which transforms need to be applied when
struct GetParsCache {
  std::unordered_set<std::string> postmap_param_set;
  std::vector<TransformSpec>      postmap_specs;
  std::vector<TransformSpec>      premap_specs;      // empty if no trend
  std::vector<TransformSpec>      pretransform_specs; // empty if no trend
};

double c_log_likelihood_race_new_path(ParamTable& pt,
                                      const RaceModelSetup& setup,
                                      const NumericVector& rts,
                                      const LogicalVector& winner,
                                      const std::vector<int>& is_ok,
                                      const std::vector<int>& idx_win,
                                      const std::vector<int>& idx_los,
                                      const IntegerVector& expand,
                                      double min_ll,
                                      int n_acc,
                                      NumericVector& raw,      // length n_trials
                                      NumericVector& ll_out,   // length idx_win.size()
                                      RaceScratch& scratch,
                                      race_fast_fn dfun_fill,  // pdf
                                      race_fast_fn pfun_fill);


// [[Rcpp::export]]
Rcpp::NumericMatrix do_transform(Rcpp::NumericMatrix pars, Rcpp::List transform) {
  // Build the specs for these parameters
  std::vector<TransformSpec> specs = make_transform_specs(pars, transform);
  // Apply transformation in place and return
  return c_do_transform(pars, specs);
}



NumericMatrix c_map_p(NumericVector p_vector,
                      CharacterVector p_types,
                      List designs,
                      int n_trials,
                      DataFrame data,
                      List trend,
                      const std::vector<TransformSpec>& full_specs,
                      bool drop_trend_pars=true) {

  // Extract information about trends
  const bool has_trend = (trend.length() > 0);
  bool premap = false;
  bool pretransform = false;
  CharacterVector trend_names;
  if (has_trend) {
    trend_names = trend.names();
    for (int i = 0; i < trend.size(); ++i) {
      List cur = trend[i];
      std::string ph = Rcpp::as<std::string>(cur["phase"]);
      if (ph == "premap") premap = true;
      if (ph == "pretransform") pretransform = true;
    }
  }

  const int n_params = p_types.size();
  NumericMatrix pars(n_trials, n_params);
  colnames(pars) = p_types;

  // Prepare trend parameter columns when needed
  NumericMatrix trend_pars;
  LogicalVector trend_index(n_params, FALSE);
  CharacterVector trend_pnames;
  if (has_trend && (premap || pretransform)) {
    // Fill in trend columns first so that they can be used in premapped trend
    // This function also applies transformations to the trend parameters
    // to ensure real-lines support.
    // The pre-transform trends are also included here, they are used after map_p
    // but they need to be transformed already (before the other trend parameters)
    // are transformed.
    trend_pars = build_trend_columns_from_design(p_vector, p_types, designs, n_trials, trend, full_specs);
    trend_pnames = colnames(trend_pars);
    trend_index = contains_multiple(p_types, trend_pnames);
  }

  // Map non-trend parameters from designs, applying premap trends if requested
  for (int i = 0; i < n_params; i++) {
    if (trend_index[i] == TRUE) continue; // skip trend parameters here
    NumericMatrix cur_design = designs[i];
    CharacterVector cur_names = colnames(cur_design);
    for (int j = 0; j < cur_design.ncol(); j++) {
      String cur_name(cur_names[j]);
      NumericVector p_mult_design(n_trials, p_vector[cur_name]);
      if (has_trend && premap) {
        p_mult_design = apply_premap_trends(data, trend, trend_names, cur_name, p_mult_design, trend_pars, p_vector);
      }
      p_mult_design = p_mult_design * cur_design(_, j);
      LogicalVector bad = is_na(p_mult_design) | is_nan(p_mult_design);
      p_mult_design[bad] = 0;
      pars(_, i) = pars(_, i) + p_mult_design;
    }
  }

  // If using pretransform trends, copy the pre-transformed trend cols into pars by name
  if (has_trend && pretransform) {
    // Only fill columns for pretransform entries
    CharacterVector tf_names = collect_trend_param_names_phase(trend, "pretransform");
    NumericMatrix trend_pars_tf = (tf_names.size() > 0) ? submat_rcpp_col_by_names(trend_pars, tf_names) : NumericMatrix(n_trials, 0);
    fill_trend_columns_for_pretransform(pars, p_types, trend_pars_tf);
  }

  // If premap, trend parameter columns are not part of the final matrix
  if (has_trend && premap && drop_trend_pars) {
    CharacterVector names_premap = collect_trend_param_names_phase(trend, "premap");
    if (names_premap.size() > 0) {
      pars = submat_rcpp_col(pars, !contains_multiple(p_types, names_premap));
    }
  }
  return pars;
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, const std::vector<PreTransformSpec>& p_specs,
                              CharacterVector p_types, List designs, int n_trials, DataFrame data, List trend,
                              const std::vector<TransformSpec>& full_specs,
                              bool drop_trend_pars=true){
  const bool has_trend = (trend.length() > 0);
  bool pretransform = false;
  bool posttransform = false;
  if (has_trend) {
    for (int i = 0; i < trend.size(); ++i) {
      List cur = trend[i];
      std::string ph = Rcpp::as<std::string>(cur["phase"]);
      if (ph == "pretransform") pretransform = true;
      if (ph == "posttransform") posttransform = true;
    }
  }
  NumericVector p_vector_updtd(clone(p_vector));
  p_vector_updtd = c_do_pre_transform(p_vector_updtd, p_specs);
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, data, trend, full_specs, drop_trend_pars);
  // // Check if pretransform trend applies
  if(pretransform){
    pars = prep_trend_phase(data, trend, pars, "pretransform");
  }
  std::vector<TransformSpec> t_specs = make_transform_specs_from_full(pars, p_types, full_specs);
  pars = c_do_transform(pars, t_specs);
  // Check if posttransform trend applies
  if(posttransform){
    // Build trend parameter columns once (transformed) and pass override
    NumericMatrix trend_pars_all = build_trend_columns_from_design(p_vector_updtd, p_types, designs, n_trials, trend, full_specs);
    CharacterVector names_post = collect_trend_param_names_phase(trend, "posttransform");
    NumericMatrix trend_pars_post = (names_post.size() > 0) ? submat_rcpp_col_by_names(trend_pars_all, names_post) : NumericMatrix(n_trials, 0);
    pars = prep_trend_phase_with_pars(data, trend, pars, "posttransform", trend_pars_post);
  }
  // ok is calculated afterwards and Ttransform applied in the function
  return(pars);
}


NumericMatrix get_pars_matrix_oo(ParamTable& param_table,
                                 const Rcpp::List& designs,
                                 TrendRuntime* trend_runtime,
                                 const GetParsCache& cache,
                                 // const std::vector<TransformSpec>& full_specs,
                                 const Rcpp::CharacterVector& keep_names,
                                 bool return_empty_matrix = false,
                                 bool return_covariate_matrix = false,
                                 bool return_all_pars = false,
                                 const std::vector<int>& kernel_output_codes = std::vector<int>{1}) {
  // Reset kernels if needed
  if (trend_runtime) {
    trend_runtime->reset_all_kernels();
  }

  // keep track of which parameters should be mapped & transformed next
  const int n_designs = designs.size();
  LogicalVector map_next(n_designs, false); // which designs we map in the *next* mapping step
  std::unordered_set<std::string> transform_next; // which params we transform in the *next* transform step

  // Trend parameter sets (may be empty if no TrendRuntime)
  // std::unordered_set<std::string> empty_set;
  // const auto& premap_set       = (trend_runtime ? trend_runtime->premap_trend_params()       : empty_set);
  // const auto& pretransform_set = (trend_runtime ? trend_runtime->pretransform_trend_params() : empty_set);

  // 3) Premap trends: MAP premap trend parameters, TRANSFORM them, RUN kernels+bases
  if (trend_runtime && trend_runtime->has_premap()) {
    // MAP: only designs whose outputs are premap trend parameters
    map_next = trend_runtime->premap_design_mask(designs);
    param_table.map_from_designs(designs, map_next);

    // TRANSFORM: only premap trend parameters
    if (!cache.premap_specs.empty()) {
      c_do_transform_pt(param_table, cache.premap_specs);
    }

    // RUN: apply all premap trend operations
    std::size_t n_ops = trend_runtime->premap_ops.size();
    for (std::size_t i = 0; i < n_ops; ++i) {
      TrendOpRuntime& op = trend_runtime->premap_ops[i];
      trend_runtime->apply_base_for_op(op, param_table);
    }
  }

  // 4) Map designs for remaining parameters. Invert map_next
  for (int i = 0; i < n_designs; ++i) {
    bool is_premap = (trend_runtime && trend_runtime->has_premap()) ? map_next[i] : false;
    map_next[i] = !is_premap;
  }
  param_table.map_from_designs(designs, map_next);

  // 5) Pretransform trends: TRANSFORM pretransform trend parameters, RUN kernels+bases
  if (trend_runtime && trend_runtime->has_pretransform()) {
    // Transform
    if (!cache.pretransform_specs.empty()) {
      c_do_transform_pt(param_table, cache.pretransform_specs);
    }

    // Trend
    std::size_t n_ops = trend_runtime->pretransform_ops.size();
    for (std::size_t i = 0; i < n_ops; ++i) {
      TrendOpRuntime& op = trend_runtime->pretransform_ops[i];
      trend_runtime->apply_base_for_op(op, param_table);
    }
  }

  // 6) Transforms for all parameters excluding the trend pars used so far.
  c_do_transform_pt(param_table, cache.postmap_specs);

  // 7) Posttransform trends.
  if (trend_runtime && trend_runtime->has_posttransform()) {
    std::size_t n_ops = trend_runtime->posttransform_ops.size();
    for (std::size_t i = 0; i < n_ops; ++i) {
      TrendOpRuntime& op = trend_runtime->posttransform_ops[i];
      trend_runtime->apply_base_for_op(op, param_table);
    }
  }

  // 8) Kernel outputs, if requested
  if (return_covariate_matrix) {
    if (!trend_runtime) {
      Rcpp::stop("return_kernel_matrix/return_covariate_matrix requested but no trend was provided");
    }

    std::vector<int> codes = kernel_output_codes;
    if (codes.empty()) {
      // reasonable default: main trajectory
      codes.push_back(1);
    }

    return trend_runtime->all_kernel_outputs(param_table, codes);
  }

  if(return_all_pars) {
    return param_table.materialize();
  }

  if(return_empty_matrix) {
    NumericMatrix tmp;
    return tmp;
  }

  // 9) Materialize only requested parameters
  return param_table.materialize_by_param_names(keep_names);
}


double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                            const int n_trials, IntegerVector expand,
                            double min_ll, std::vector<int> is_ok){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  NumericVector lls(n_trials);
  lls = d_DDM_Wien(rts, R, pars, is_ok);

  // lls_exp = c_expand(lls, expand); // decompress
  // // lls_exp = lls;
  // lls_exp[is_na(lls_exp)] = min_ll;
  // lls_exp[is_infinite(lls_exp)] = min_ll;
  // lls_exp[lls_exp < min_ll] = min_ll;
  // return(sum(lls_exp));
  // More SIMD-friendly == faster
  // decompress

  const double* lls_ptr    = lls.begin();
  const int*    expand_ptr = expand.begin();

  double sum_ll = 0.0;

  // expand is 1-based, so subtract 1
#pragma omp simd reduction(+:sum_ll)
  for (int i = 0; i < n_out; ++i) {
    int idx = expand_ptr[i] - 1;
    double v = lls_ptr[idx];

    if (!std::isfinite(v) || v < min_ll) {
      v = min_ll;
    }
    sum_ll += v;
  }

  return sum_ll;
}


double c_log_likelihood_race(NumericMatrix pars, DataFrame data,
                             NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                             NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                             const int n_trials, LogicalVector winner, IntegerVector expand,
                             double min_ll, LogicalVector is_ok){
  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  CharacterVector R = data["R"];
  NumericVector lR = data["lR"];
  NumericVector lds_exp(n_out);
  const int n_acc = unique(lR).length();
  if(sum(contains(data.names(), "RACE")) == 1){
    NumericVector NACC = data["RACE"];
    CharacterVector vals_NACC = NACC.attr("levels");
    for(int x = 0; x < pars.nrow(); x++){
      // subtract 1 because R is 1 coded
      if(lR[x] > atoi(vals_NACC[NACC[x]-1])){
        pars(x,0) = NA_REAL;
      }
    }
  }

  // Precompute winner / loser indices once
  int* win_flag = LOGICAL(winner);
  std::vector<int> idx_win;
  std::vector<int> idx_los;
  idx_win.reserve(n_trials);
  idx_los.reserve(n_trials);

  for (int i = 0; i < n_trials; ++i) {
    if (win_flag[i])
      idx_win.push_back(i);
    else
      idx_los.push_back(i);
  }

  // 1) Winner densities -> log -> into lds at winner positions
  NumericVector win = log(dfun(rts, pars, winner, std::exp(min_ll), is_ok));

  if ((int)win.size() != (int)idx_win.size()) {
    Rcpp::stop("c_log_likelihood_race: dfun() length mismatch with winner mask");
  }

  double* lds_ptr = lds.begin();
  double* win_ptr = win.begin();

  for (int k = 0; k < (int)idx_win.size(); ++k) {
    lds_ptr[ idx_win[k] ] = win_ptr[k];
  }

  // 2) Loser side (if more than one accumulator)
  if (n_acc > 1) {
    NumericVector loss = log(1 - pfun(rts, pars, !winner, std::exp(min_ll), is_ok));

    if ((int)loss.size() != (int)idx_los.size()) {
      Rcpp::stop("c_log_likelihood_race: pfun() length mismatch with !winner mask");
    }

    double* loss_ptr = loss.begin();
    const double sentinel = std::log(1.0 - std::exp(min_ll));  // log(1 - exp(min_ll))

    for (int k = 0; k < (int)idx_los.size(); ++k) {
      double v = loss_ptr[k];
      // original behaviour:
      // loss[is_na(loss)] = min_ll;
      // loss[loss == log(1 - exp(min_ll))] = min_ll;
      if (!std::isfinite(v) || v == sentinel) {
        v = min_ll;
      }
      loss_ptr[k] = v;
      lds_ptr[ idx_los[k] ] = v;
    }
  }

  // 3) Replace any remaining NA/Inf in lds with min_ll (from win side)
  for (int i = 0; i < n_trials; ++i) {
    double v = lds_ptr[i];
    if (!std::isfinite(v)) {
      lds_ptr[i] = min_ll;
    }
  }
  // 4) Combine across accumulators and expand

  if (n_acc > 1) {
    // We expect one winner per "trial group"; n_winner == number of groups here.
    const int n_winners = static_cast<int>(idx_win.size());
    const int n_losers  = static_cast<int>(idx_los.size());

    if (n_losers != n_winners * (n_acc - 1)) {
      Rcpp::stop("c_log_likelihood_race: n_losers != n_winners * (n_acc - 1)");
    }

    // Build ll_out (one per winner) and lds_los (all loser entries)
    NumericVector ll_out(n_winners);
    NumericVector lds_los(n_losers);
    double* ll_ptr  = ll_out.begin();
    double* los_ptr = lds_los.begin();

    // Fill winners
    for (int t = 0; t < n_winners; ++t) {
      ll_ptr[t] = lds_ptr[ idx_win[t] ];
    }

    // Fill losers in order
    for (int k = 0; k < n_losers; ++k) {
      los_ptr[k] = lds_ptr[ idx_los[k] ];
    }

    // Add loser contributions
    if (n_acc == 2) {
      // One loser per winner; simple pairwise add
      for (int t = 0; t < n_winners; ++t) {
        ll_ptr[t] += los_ptr[t];
      }
    } else {
      const int stride = n_acc - 1;
      for (int t = 0; t < n_winners; ++t) {
        double s = 0.0;
        const int base = t * stride;
        for (int k = 0; k < stride; ++k) {
          s += los_ptr[base + k];
        }
        ll_ptr[t] += s;
      }
    }

    // Decompress and clamp+sum
    NumericVector ll_exp = c_expand(ll_out, expand);
    double* x    = ll_exp.begin();
    const int m  = ll_exp.size();
    double sum_ll = 0.0;

#pragma omp simd reduction(+:sum_ll)
    for (int i = 0; i < m; ++i) {
      double v = x[i];
      if (!std::isfinite(v) || v < min_ll) {
        v = min_ll;
      }
      // x[i] = v;
      sum_ll += v;
    }
    return sum_ll;

  } else {
    // Single accumulator: just expand lds and clamp+sum
    NumericVector lds_exp = c_expand(lds, expand);
    double* x    = lds_exp.begin();
    const int m  = lds_exp.size();
    double sum_ll = 0.0;

#pragma omp simd reduction(+:sum_ll)
    for (int i = 0; i < m; ++i) {
      double v = x[i];
      if (!std::isfinite(v) || v < min_ll) {
        v = min_ll;
      }
      // x[i] = v;
      sum_ll += v;
    }
    return sum_ll;
  }
}



// //#// [[Rcpp::export]]
// NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
//                       List designs, String type, List bounds, List transforms, List pretransforms,
//                       CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue){ //,
// //            bool debug_first_particle = false){
//   const int n_particles = p_matrix.nrow();
//   const int n_trials = data.nrow();
//   NumericVector lls(n_particles);
//   NumericVector p_vector(p_matrix.ncol());
//   CharacterVector p_names = colnames(p_matrix);
//   p_vector.names() = p_names;
//   NumericMatrix pars;
//   LogicalVector is_ok(n_trials);
//
//   // Once (outside the main loop over particles):
//   NumericMatrix minmax = bounds["minmax"];
//   CharacterVector mm_names = colnames(minmax);
//   std::vector<PreTransformSpec> p_specs;
//   std::vector<BoundSpec> bound_specs;
//   std::vector<TransformSpec> full_t_specs; // precomputed transform specs for p_types
//
//   /// build trend plan, runtime (only needed for new pathway), and paramtable
//   NumericVector dummy_p = clone(p_vector);
//   ParamTable pt_template = ParamTable::from_p_vector_and_designs(dummy_p, designs, data.nrow());
//   std::vector<TransformSpec> full_specs_new;
//
//   Rcpp::List trend_list;
//
//   if (!trend.isNull()) {
//     trend_list = Rcpp::List(trend);   // safe: wraps the underlying SEXP
//   } else {
//     trend_list = Rcpp::List::create();  // or Rcpp::List(); empty list
//   }
//   if(type == "DDM"){
//     IntegerVector expand = data.attr("expand");
//     for(int i = 0; i < n_particles; i++){
//       p_vector = p_matrix(i, _);
//       if(i == 0){
//         p_specs = make_pretransform_specs(p_vector, pretransforms);
//         // Precompute transform specs for all p_types using a one-time dummy
//         NumericMatrix dummy(1, p_types.size());
//         colnames(dummy) = p_types;
//         full_t_specs = make_transform_specs(dummy, transforms);
//       }
//         pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend_list, full_t_specs);
//       // Precompute specs
//       if (i == 0) {                            // first particle only, just to get colnames
//         bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
//       }
//       is_ok = c_do_bound(pars, bound_specs);
//       lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
//     }
//   } else if(type == "MRI" || type == "MRI_AR1"){
//     int n_pars = p_types.length();
//     NumericVector y = extract_y(data);
//     for(int i = 0; i < n_particles; i++){
//       p_vector = p_matrix(i, _);
//       if(i == 0){
//         p_specs = make_pretransform_specs(p_vector, pretransforms);
//         // Precompute transform specs for all p_types using a one-time dummy
//         NumericMatrix dummy(1, p_types.size());
//         colnames(dummy) = p_types;
//         full_t_specs = make_transform_specs(dummy, transforms);
//       }
//         pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend_list, full_t_specs);
//       // Precompute specs
//       if (i == 0) {                            // first particle only, just to get colnames
//         bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
//       }
//       is_ok = c_do_bound(pars, bound_specs);
//       if(type == "MRI"){
//         lls[i] = c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
//       } else{
//         lls[i] = c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
//       }
//     }
//   } else{
//     IntegerVector expand = data.attr("expand");
//     LogicalVector winner = data["winner"];
//     // Love me some good old ugly but fast c++ pointers
//     NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
//     NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
//     if(type == "LBA"){
//       dfun = dlba_c;
//       pfun = plba_c;
//     } else if(type == "RDM"){
//       dfun = drdm_c;
//       pfun = prdm_c;
//     } else {
//       dfun = dlnr_c;
//       pfun = plnr_c;
//     }
//     NumericVector lR = data["lR"];
//     int n_lR = unique(lR).length();
//     for (int i = 0; i < n_particles; ++i) {
//       p_vector = p_matrix(i, _);
//       if(i == 0){
//         p_specs = make_pretransform_specs(p_vector, pretransforms);
//         // Precompute transform specs for all p_types using a one-time dummy
//         NumericMatrix dummy(1, p_types.size());
//         colnames(dummy) = p_types;
//         full_t_specs = make_transform_specs(dummy, transforms);
//       }
//
//       pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend_list, full_t_specs);
//       if (i == 0) {                            // first particle only, just to get colnames
//         bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
//       }
//       is_ok = c_do_bound(pars, bound_specs);
//       is_ok = lr_all(is_ok, n_lR);
//       lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
//     }
//   }
//   return(lls);
// }



// [[Rcpp::export]]
NumericVector calc_ll_oo(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                         List designs, String type, List bounds, List transforms, List pretransforms,
                         CharacterVector p_types, double min_ll,
                         Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                         bool use_pt_pipeline = true)
{
  const int n_particles = particle_matrix.nrow();
  const int n_trials    = data.nrow();

  NumericMatrix  pars;
  // LogicalVector  is_ok(n_trials);
  NumericVector  lls(n_particles);
  std::vector<int> is_ok(n_trials, 1);

  // --- Bounds ---
  NumericMatrix   minmax   = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs;

  // 1. Pre-transform particle matrix
  std::vector<TransformSpec> t_specs = make_transform_specs(particle_matrix, pretransforms);
  particle_matrix = c_do_transform(particle_matrix, t_specs);

  // 2. Append constants (once)
  bool has_constants = !(constants.size() == 1 && Rcpp::NumericVector::is_na(constants[0]));
  if (has_constants) {
    particle_matrix = add_constants_columns(particle_matrix, constants);
  }

  // 3. Build ParamTable template from first particle
  NumericVector p_vector = particle_matrix(0, _);
  p_vector.attr("names") = colnames(particle_matrix);
  ParamTable param_table_template = ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials);

  // 4. Transform specs for ParamTable
  std::vector<TransformSpec> transform_specs =
    make_transform_specs_for_paramtable(param_table_template, transforms);

  // 5. Trend objects and keep_names
  std::unique_ptr<TrendPlan>    trend_plan;
  std::unique_ptr<TrendRuntime> trend_runtime;
  Rcpp::CharacterVector keep_names;

  if (!trend.isNull()) {
    trend_plan.reset(new TrendPlan(trend, data));
    trend_runtime.reset(new TrendRuntime(*trend_plan));
    trend_runtime->bind_all_ops_to_paramtable(param_table_template);

    Rcpp::CharacterVector dnames = designs.names();
    const auto& trend_params = trend_runtime->all_trend_params();
    keep_names = names_excluding(dnames, { &trend_params });
  } else {
    keep_names = designs.names();
  }

  TrendRuntime* trend_runtime_ptr = trend_runtime ? trend_runtime.get() : nullptr;

  // 6. Fast column-index lookup: particle matrix column → ParamTable base index
  Rcpp::CharacterVector pm_names = colnames(particle_matrix);
  std::vector<int> pm_col_to_base_idx(pm_names.size(), -1);
  for (int j = 0; j < pm_names.size(); ++j) {
    std::string nm = Rcpp::as<std::string>(pm_names[j]);
    auto it = param_table_template.name_to_base_idx.find(nm);
    if (it != param_table_template.name_to_base_idx.end()) {
      pm_col_to_base_idx[j] = it->second;
    }
  }

  // Cache which parameters need to be transformed at which point
  static const std::unordered_set<std::string> empty_set;  // static: constructed once, never changes
  GetParsCache gp_cache;
  {
    const auto& premap_set       = trend_runtime_ptr ? trend_runtime_ptr->premap_trend_params()       : empty_set;
    const auto& pretransform_set = trend_runtime_ptr ? trend_runtime_ptr->pretransform_trend_params() : empty_set;

    gp_cache.postmap_param_set = param_names_excluding(param_table_template, { &premap_set, &pretransform_set });
    gp_cache.postmap_specs     = filter_specs_by_param_set(param_table_template, transform_specs, gp_cache.postmap_param_set);

    if (trend_runtime_ptr && trend_runtime_ptr->has_premap()) {
      gp_cache.premap_specs = filter_specs_by_param_set(param_table_template, transform_specs, premap_set);
    }
    if (trend_runtime_ptr && trend_runtime_ptr->has_pretransform()) {
      gp_cache.pretransform_specs = filter_specs_by_param_set(param_table_template, transform_specs, pretransform_set);
    }
  }

  // -----------------------------------------------------------------------
  // DDM
  // -----------------------------------------------------------------------
  if (type == "DDM") {
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) param_table_template.fill_from_particle_row(particle_matrix, i, pm_col_to_base_idx);
      pars = get_pars_matrix_oo(param_table_template, designs, trend_runtime_ptr,
                                gp_cache, keep_names);
      if (i == 0) bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      c_do_bound_pt(param_table_template, bound_specs, is_ok);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
    }

  // -----------------------------------------------------------------------
  // MRI / MRI_AR1
  // -----------------------------------------------------------------------
  } else if (type == "MRI" || type == "MRI_AR1") {
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) param_table_template.fill_from_particle_row(particle_matrix, i, pm_col_to_base_idx);
      pars = get_pars_matrix_oo(param_table_template, designs, trend_runtime_ptr,
                                gp_cache, keep_names);
      if (i == 0) bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      // is_ok  = c_do_bound_pt(param_table_template, bound_specs);
      c_do_bound_pt(param_table_template, bound_specs, is_ok);
      lls[i] = (type == "MRI")
        ? c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll)
          : c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
    }

  // -----------------------------------------------------------------------
  // Race models (RDM, LBA, LNR)
  // -----------------------------------------------------------------------
  } else {
    NumericVector lR     = data["lR"];
    IntegerVector expand = data.attr("expand");
    const int n_acc      = unique(lR).length();
    NumericVector rts    = data["rt"];
    LogicalVector winner = data["winner"];

    // Precompute winner/loser index lists (once, outside particle loop)
    std::vector<int> idx_win, idx_los;
    idx_win.reserve(n_trials);
    idx_los.reserve(n_trials);
    int* win_flag = LOGICAL(winner);
    for (int i = 0; i < n_trials; ++i) {
      (win_flag[i] ? idx_win : idx_los).push_back(i);
    }
    const int n_winners = (int)idx_win.size();

    // Scratch buffers (reused across particles)
    NumericVector raw(n_trials);
    NumericVector ll_out(n_winners);

    RaceScratch scratch;
    scratch.reserve(std::max((int)idx_win.size(), (int)idx_los.size()));

    // Race model setup
    RaceModelSetup setup = make_race_setup(type, param_table_template);

    // Variable-accumulator (RACE column) support — hoist out of particle loop
    const bool has_race_col = (sum(contains(data.names(), "RACE")) == 1);
    NumericVector   NACC;
    CharacterVector vals_NACC;
    if (has_race_col) {
      NACC      = data["RACE"];
      vals_NACC = NACC.attr("levels");
    }

    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) param_table_template.fill_from_particle_row(particle_matrix, i, pm_col_to_base_idx);
      get_pars_matrix_oo(param_table_template, designs, trend_runtime_ptr,
                         gp_cache, keep_names,
                         /*return_empty_matrix=*/true);

      if (i == 0) bound_specs = make_bound_specs_pt(minmax, mm_names, param_table_template, bounds);
      c_do_bound_pt(param_table_template, bound_specs, is_ok);
      // is_ok = c_do_bound_pt(param_table_template, bound_specs);
      lr_all(is_ok, n_acc);

      if (has_race_col) {
        NumericMatrix& base = param_table_template.base;
        for (int x = 0; x < base.nrow(); ++x) {
          if (lR[x] > atoi(vals_NACC[NACC[x] - 1])) {
            base(x, setup.col_na_marker) = NA_REAL;
          }
        }
      }

      lls[i] = c_log_likelihood_race_new_path(
        param_table_template, setup,
        rts, winner, is_ok,
        idx_win, idx_los, expand,
        min_ll, n_acc, raw, ll_out,
        scratch,
        setup.fill_pdf, setup.fill_cdf);

    }
  }

  return lls;
}


// [[Rcpp::export]]
NumericMatrix get_pars_c_wrapper(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
                                 List designs, List bounds, List transforms, List pretransforms,
                                 CharacterVector p_types,
                                 Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                 bool return_kernel_matrix = false,
                                 bool drop_trend_pars = true)
{

  if (return_kernel_matrix && trend.isNull()) {
    stop("return_kernel_matrix=TRUE requires a non-NULL 'trend' specification");
  }

  const int n_trials = data.nrow();

  // take first row as in your original code
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  p_vector.names() = p_names;
  p_vector = p_matrix(0, _);

  NumericMatrix pars;

  // Pretransform specs and transform specs (once)
  std::vector<PreTransformSpec> p_specs;
  std::vector<TransformSpec> full_t_specs;

  p_specs = make_pretransform_specs(p_vector, pretransforms);

  NumericMatrix dummy(1, p_types.size());
  colnames(dummy) = p_types;
  full_t_specs = make_transform_specs(dummy, transforms);


  // Old pathway: needs a List, not Nullable
  Rcpp::List trend_list;
  if (!trend.isNull()) {
    trend_list = Rcpp::List(trend);   // wraps underlying SEXP
  } else {
    trend_list = Rcpp::List::create();  // empty list => no trend
  }

  pars = get_pars_matrix(p_vector, constants,
                         p_specs, p_types,
                         designs, n_trials, data,
                         trend_list, full_t_specs, drop_trend_pars);

  return pars;
}

// [[Rcpp::export]]
NumericMatrix get_pars_c_wrapper_oo(NumericMatrix particle_matrix,
                                    DataFrame data,
                                    NumericVector constants,
                                    List designs,
                                    List bounds,
                                    List transforms,
                                    List pretransforms,
                                    Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                    bool return_kernel_matrix = false,
                                    bool return_all_pars = false,
                                    IntegerVector kernel_output_codes = 1)
{
  const int n_trials = data.nrow();

  if (Rf_isNull(colnames(particle_matrix))) {
    stop("p_matrix must have column names for pretransforms/transform specs");
  }

  // 1. Pre-transform p_matrix - using old pipeline
  std::vector<TransformSpec> t_specs = make_transform_specs(particle_matrix, pretransforms);
  particle_matrix = c_do_transform(particle_matrix, t_specs);

  // 2. Add constants. Makes a copy - not ideal but hopefully just a minor little bit of overhead
  // constants: NumericVector (may be a single NA meaning "no constants")
  bool has_constants = !(constants.size() == 1 &&
                         Rcpp::NumericVector::is_na(constants[0]));
  if (has_constants) {
    particle_matrix = add_constants_columns(particle_matrix, constants);
  }

  // 3. Start building objects needed for get_pars_matrix.
  // 3.1 Param_table
  NumericVector p_vector = particle_matrix(0, _);
  p_vector.attr("names") = colnames(particle_matrix);
  ParamTable param_table_template = ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials);

  // 3.2 Transform specs
  std::vector<TransformSpec> transform_specs;
  transform_specs = make_transform_specs_for_paramtable(param_table_template, transforms);

  // 3.3 Trend objects and return_param_names (which parameters to return)
  std::unique_ptr<TrendPlan>    trend_plan;
  std::unique_ptr<TrendRuntime> trend_runtime;
  Rcpp::CharacterVector return_param_names;

  if (!trend.isNull()) {
    // Build TrendPlan/TrendRuntime only if trend is provided
    trend_plan.reset(new TrendPlan(trend, data));
    trend_runtime.reset(new TrendRuntime(*trend_plan));

    // Bind ops to the fixed ParamTable layout once
    trend_runtime->bind_all_ops_to_paramtable(param_table_template);

    // return_param_names = design names minus trend params
    Rcpp::CharacterVector dnames = designs.names();
    const auto& trend_params = trend_runtime->all_trend_params();
    return_param_names = names_excluding(dnames, { &trend_params });
  } else {
    // No trend: keep all design parameters
    Rcpp::CharacterVector dnames = designs.names();
    if (Rf_isNull(dnames)) {
      // Either keep nothing, or better, derive from ParamTable
      // e.g. keep all parameters that come from designs:
      // keep_names = param_table_template.design_param_names();
      return_param_names = Rcpp::CharacterVector(0);  // at least STRSXP, not NULL
    } else {
      return_param_names = dnames;
    }
    // trend_runtime stays null
  }

  TrendRuntime* trend_runtime_ptr = trend_runtime ? trend_runtime.get() : nullptr;



  // Convert IntegerVector -> std::vector<int>
  std::vector<int> kernel_codes;
  kernel_codes.reserve(kernel_output_codes.size());
  for (int i = 0; i < kernel_output_codes.size(); ++i) {
    kernel_codes.push_back(kernel_output_codes[i]);
  }

  // get pars cache
  static const std::unordered_set<std::string> empty_set;  // static: constructed once, never changes
  GetParsCache gp_cache;
  {
    const auto& premap_set       = trend_runtime_ptr ? trend_runtime_ptr->premap_trend_params()       : empty_set;
    const auto& pretransform_set = trend_runtime_ptr ? trend_runtime_ptr->pretransform_trend_params() : empty_set;

    gp_cache.postmap_param_set = param_names_excluding(param_table_template, { &premap_set, &pretransform_set });
    gp_cache.postmap_specs     = filter_specs_by_param_set(param_table_template, transform_specs, gp_cache.postmap_param_set);

    if (trend_runtime_ptr && trend_runtime_ptr->has_premap()) {
      gp_cache.premap_specs = filter_specs_by_param_set(param_table_template, transform_specs, premap_set);
    }
    if (trend_runtime_ptr && trend_runtime_ptr->has_pretransform()) {
      gp_cache.pretransform_specs = filter_specs_by_param_set(param_table_template, transform_specs, pretransform_set);
    }
  }

  NumericMatrix pars = get_pars_matrix_oo(param_table_template,
                                          designs,
                                          trend_runtime_ptr,
                                          gp_cache,
                                          return_param_names,
                                          false, // return_empty_matrix - set to false
                                          return_kernel_matrix,
                                          return_all_pars,
                                          kernel_codes);
  return pars;
}


/// New pathway
// ---------------------------------------------------------------------------
// c_log_likelihood_race_new_path
// ---------------------------------------------------------------------------

double c_log_likelihood_race_new_path(ParamTable& pt,
                                      const RaceModelSetup& setup,
                                      const NumericVector& rts,
                                      const LogicalVector& winner,
                                      const std::vector<int>& is_ok,
                                      const std::vector<int>& idx_win,
                                      const std::vector<int>& idx_los,
                                      const IntegerVector& expand,
                                      double min_ll,
                                      int n_acc,
                                      NumericVector& raw,
                                      NumericVector& ll_out,
                                      RaceScratch& scratch,
                                      race_fast_fn dfun_fill,
                                      race_fast_fn pfun_fill)
{
  const int n_winners = (int)idx_win.size();
  const int n_losers  = (int)idx_los.size();

  double* raw_ptr = raw.begin();
  double* ll_ptr  = ll_out.begin();
  const int* ok_ptr = is_ok.data();

  // 1) Fill log(pdf) for winners and log(1-cdf) for losers into raw.
  //    fill_both stores pdf / (1-cdf); vec_log transforms the whole array
  //    in one vectorised pass (vvlog on Apple, libmvec on Linux/x86).
  //    Invalid inputs (<=0, nan) produce -inf or nan, which the clamp below
  //    catches — no per-element branching needed.
  //
  //   // setup.fill_both() refers to gather-scatter implementations.
  //   // on linux/x86, this is significantly faster. macOS/arm64 doesn't care

  setup.fill_both(rts, pt, setup.spec, idx_win, idx_los, raw_ptr, scratch);
  vec_log(raw_ptr, raw.size());  // bulk log over entire raw buffer

  // 2) Per-trial log-likelihood into ll_out.
  //    raw now contains log(pdf) at winner indices, log(1-cdf) at loser indices.
  //    Clamp to min_ll using !(v > min_ll) which catches -inf and nan.
  auto clamp = [min_ll](double v) {
    return (v > min_ll) ? v : min_ll;
  };

  if (n_acc == 1) {
    for (int t = 0; t < n_winners; ++t) {
      const int i_win = idx_win[t];
      ll_ptr[t] = ok_ptr[i_win] ? clamp(raw_ptr[i_win]) : min_ll;
    }
  } else {
    const int stride = n_acc - 1;
    if (n_losers != n_winners * stride) {
      Rcpp::stop("c_log_likelihood_race_new_path: n_losers != n_winners * (n_acc - 1)");
    }

    for (int t = 0; t < n_winners; ++t) {
      const int i_win = idx_win[t];

      //
      if (!ok_ptr[i_win]) {
        ll_ptr[t] = min_ll;
        continue;
      }

      double ll = clamp(raw_ptr[i_win]);

      const int base = t * stride;
      for (int k = 0; k < stride; ++k) {
        ll += clamp(raw_ptr[idx_los[base + k]]);  // no ok check needed, ok for winners and losers are guaranteed to match by lr_all()
      }

      ll_ptr[t] = ll;
    }
  }

  // 3) Expand and sum
  const int  m       = expand.size();
  const int* exp_ptr = expand.begin();
  double sum_ll = 0.0;

#pragma omp simd reduction(+:sum_ll)
  for (int i = 0; i < m; ++i) {
    sum_ll += clamp(ll_ptr[exp_ptr[i] - 1]);
  }

//   // 1) Fill pdfs (winners) and cdfs (losers)
//   // dfun_fill(rts, pt, setup.spec, winner, raw_ptr);
//   // if (n_acc > 1) {
//   //   pfun_fill(rts, pt, setup.spec, winner, raw_ptr);
//   // }
//
//   // setup.fill_both() refers to gather-scatter implementations.
//   // on linux/x86, this is significantly faster. macOS/arm64 doesn't care
//   setup.fill_both(rts, pt, setup.spec, idx_win, idx_los, raw_ptr, scratch);
//
//   // 2) Per-trial log-likelihood into ll_out
//   // safe_log: returns log(x) clamped to [min_ll, inf), or min_ll if x <= 0
//   // Defined as a static helper to stay C++11 compatible (no generic lambda)
//   struct SafeLog {
//     double min_ll;
//     double operator()(double x) const {
//       if (!std::isfinite(x) || x <= 0.0) return min_ll;
//       double v = std::log(x);
//       return (std::isfinite(v) && v >= min_ll) ? v : min_ll;
//     }
//   } safe_log = { min_ll };
//
//   if (n_acc == 1) {
//     for (int t = 0; t < n_winners; ++t) {
//       const int i_win = idx_win[t];
//       ll_ptr[t] = ok_ptr[i_win] ? safe_log(raw_ptr[i_win]) : min_ll;
//     }
//   } else {
//     const int stride = n_acc - 1;
//     if (n_losers != n_winners * stride) {
//       Rcpp::stop("c_log_likelihood_race_new_path: n_losers != n_winners * (n_acc - 1)");
//     }
//
//     for (int t = 0; t < n_winners; ++t) {
//       const int i_win = idx_win[t];
//
//       if (!ok_ptr[i_win]) {
//         ll_ptr[t] = min_ll;
//         continue;
//       }
//
//       double ll = safe_log(raw_ptr[i_win]);
//
//       const int base = t * stride;
//       for (int k = 0; k < stride; ++k) {
//         const int i_los = idx_los[base + k];
//         if (!ok_ptr[i_los]) {
//           ll += min_ll;
//           continue;
//         }
//         double cdf = raw_ptr[i_los];
//         if (!std::isfinite(cdf)) {
//           ll += min_ll;
//           continue;
//         }
//         double one_minus = cdf; //1.0 - cdf;
//         ll += (one_minus <= 0.0) ? min_ll : safe_log(one_minus);
//       }
//
//       ll_ptr[t] = ll;
//     }
//   }
//
//   // 3) Expand and sum
//   const int  m       = expand.size();
//   const int* exp_ptr = expand.begin();
//   double sum_ll = 0.0;
//
// #pragma omp simd reduction(+:sum_ll)
//   for (int i = 0; i < m; ++i) {
//     double v = ll_ptr[exp_ptr[i] - 1];
//     sum_ll += (std::isfinite(v) && v >= min_ll) ? v : min_ll;
//   }

  return sum_ll;
}
// double c_log_likelihood_race_new_path(ParamTable& pt,
//                                       const void* model_spec,
//                                       const NumericVector& rts,
//                                       const LogicalVector& winner,
//                                       const LogicalVector& is_ok,
//                                       const std::vector<int>& idx_win,
//                                       const std::vector<int>& idx_los,
//                                       const IntegerVector& expand,
//                                       double min_ll,
//                                       int n_acc,
//                                       NumericVector& raw,      // length n_trials
//                                       NumericVector& ll_out,   // length idx_win.size()
//                                       race_fast_fn dfun_fill,  // pdf fill
//                                       race_fast_fn pfun_fill)  // cdf fill
// {
//   const int n_winners = (int)idx_win.size();
//   const int n_losers  = (int)idx_los.size();
//
//   double* raw_ptr = raw.begin();
//   double* ll_ptr  = ll_out.begin();
//
//   // 1) Fill pdfs/cdfs using the winner mask
//   dfun_fill(rts, pt, model_spec, winner, raw_ptr);   // winners: pdfs in raw[i]
//   if (n_acc > 1) {
//     pfun_fill(rts, pt, model_spec, winner, raw_ptr); // losers: cdfs in raw[i]
//   }
//
//   int* ok_ptr = LOGICAL(is_ok);
//
//   // 2) Per-winner loglik (winner + its losers) into ll_out
//   if (n_acc == 1) {
//     for (int t = 0; t < n_winners; ++t) {
//       const int i_win = idx_win[t];
//       double ll;
//       if (!ok_ptr[i_win]) {
//         ll = min_ll;
//       } else {
//         double pdf = raw_ptr[i_win];
//         if (!std::isfinite(pdf) || pdf <= 0.0) {
//           ll = min_ll;
//         } else {
//           ll = std::log(pdf);
//           if (!std::isfinite(ll) || ll < min_ll) ll = min_ll;
//         }
//       }
//       ll_ptr[t] = ll;
//     }
//   } else {
//     const int stride = n_acc - 1;
//     if (n_losers != n_winners * stride) {
//       Rcpp::stop("c_log_likelihood_race_new_path: n_losers != n_winners * (n_acc - 1)");
//     }
//
//     for (int t = 0; t < n_winners; ++t) {
//       const int i_win = idx_win[t];
//
//       double ll;
//       if (!ok_ptr[i_win]) {
//         ll = min_ll;
//       } else {
//         // winner part
//         double pdf = raw_ptr[i_win];
//         if (!std::isfinite(pdf) || pdf <= 0.0) {
//           ll = min_ll;
//         } else {
//           ll = std::log(pdf);
//           if (!std::isfinite(ll) || ll < min_ll) ll = min_ll;
//         }
//
//         // loser contributions
//         const int base = t * stride;
//         for (int k = 0; k < stride; ++k) {
//           const int i_los = idx_los[base + k];
//           double term;
//           if (!ok_ptr[i_los]) {
//             term = min_ll;
//           } else {
//             double cdf = raw_ptr[i_los];
//             if (!std::isfinite(cdf)) {
//               term = min_ll;
//             } else {
//               double one_minus = 1.0 - cdf;
//               if (one_minus <= 0.0) {
//                 term = min_ll;
//               } else {
//                 term = std::log(one_minus);
//                 if (!std::isfinite(term) || term < min_ll) term = min_ll;
//               }
//             }
//           }
//           ll += term;
//         }
//       }
//       ll_ptr[t] = ll;
//     }
//   }
//
//   // 3) Expand and clamp+sum
//   const int m        = expand.size();
//   const int* exp_ptr = expand.begin();
//
//   double sum_ll = 0.0;
//
// #pragma omp simd reduction(+:sum_ll)
//   for (int i = 0; i < m; ++i) {
//     int idx = exp_ptr[i] - 1;
//     double v = ll_ptr[idx];
//     if (!std::isfinite(v) || v < min_ll) {
//       v = min_ll;
//     }
//     sum_ll += v;
//   }
//
//   return sum_ll;
// }

