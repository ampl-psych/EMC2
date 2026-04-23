// Only debugging ode

// #include <Rcpp.h>
// // #include <RcppArmadillo.h>
// #include "ParamTable.h"
// #include "TrendEngine.h"   // contains TrendEngine definition
// #include "transform_utils.h"
// using namespace Rcpp;
//
// // Lightweight wrapper for R-facing API: holds both plan and runtime
// struct TrendEngine {
//   TrendPlan plan;
//   TrendRuntime runtime;
//
//   TrendEngine(Rcpp::Nullable<Rcpp::List> trend_,
//               const Rcpp::DataFrame& data_)
//     : plan(trend_, data_), runtime(plan) {}
// };
//
// // Factory: create a ParamTable from p_types and n_trials
// // [[Rcpp::export]]
// SEXP ParamTable_create_from_p_types(int n_trials,
//                                     CharacterVector p_types) {
//   ParamTable* pt = new ParamTable(ParamTable::from_p_types(n_trials, p_types));
//   XPtr<ParamTable> ptr(pt, true); // true => R will delete pt when GC'ed
//   return ptr;
// }
//
// // Factory: create from an existing matrix + names (useful for debugging)
// // [[Rcpp::export]]
// SEXP ParamTable_create_from_matrix(NumericMatrix base,
//                                    CharacterVector names) {
//   ParamTable* pt = new ParamTable(base, names);
//   XPtr<ParamTable> ptr(pt, true);
//   return ptr;
// }
//
// // Wrapper: call materialize()
// // [[Rcpp::export]]
// NumericMatrix ParamTable_materialize(SEXP pt_xptr) {
//   XPtr<ParamTable> pt(pt_xptr);
//   return pt->materialize();
// }
//
// // // Wrapper: drop columns by name
// // // [[Rcpp::export]]
// // void ParamTable_drop(SEXP pt_xptr, CharacterVector drop_names) {
// //   XPtr<ParamTable> pt(pt_xptr);
// //   pt->drop(drop_names);
// // }
//
// // Wrapper: set a column by name
// // [[Rcpp::export]]
// void ParamTable_set_column(SEXP pt_xptr,
//                            std::string name,
//                            NumericVector col) {
//   XPtr<ParamTable> pt(pt_xptr);
//   pt->set_column_by_name(name, col);
// }
//
// // Wrapper
// // [[Rcpp::export]]
// SEXP ParamTable_create_from_pvector_designs(NumericVector p_vector,
//                                             List designs,
//                                             int n_trials) {
//   ParamTable* pt = new ParamTable(
//     ParamTable::from_p_vector_and_designs(p_vector, designs, n_trials)
//   );
//   XPtr<ParamTable> ptr(pt, true);
//   return ptr;
// }
//
//
// // [[Rcpp::export]]
// void ParamTable_map_designs(SEXP pt_xptr,
//                             Rcpp::List designs,
//                             Rcpp::LogicalVector include_param) {
//   Rcpp::XPtr<ParamTable> pt(pt_xptr);
//   pt->map_from_designs(designs, include_param);
// }
//
// // // [[Rcpp::export]]
// // SEXP TrendEngine_create(List trend, DataFrame data) {
// //   // allocate on heap and wrap in XPtr so R manages lifetime
// //   TrendEngine* engine = new TrendEngine(trend, data);
// //   Rcpp::XPtr<TrendEngine> xp(engine, true); // true => delete when GC'd
// //   return xp;
// // }
// // [[Rcpp::export]]
// SEXP TrendEngine_create(SEXP trend, Rcpp::DataFrame data) {
//   Rcpp::Nullable<Rcpp::List> ntrend(trend);
//   TrendEngine* engine = new TrendEngine(ntrend, data);
//   Rcpp::XPtr<TrendEngine> xp(engine, true);
//   return xp;
// }
//
// // [[Rcpp::export]]
// LogicalVector TrendEngine_premap_mask(SEXP engine_ptr, List designs) {
//   Rcpp::XPtr<TrendEngine> engine(engine_ptr);
//   return engine->plan.premap_design_mask(designs);
// }
//
//
// // [[Rcpp::export]]
// int ParamTable_get_n_trials(SEXP pt_xptr) {
//   XPtr<ParamTable> pt(pt_xptr);
//   return pt->n_trials;   // or a getter if you made n_trials private
// }
//
//
// // [[Rcpp::export]]
// void ParamTable_do_transform_all(SEXP param_table_ptr, Rcpp::List transform) {
//   // however you're storing/accessing ParamTable; here assume by reference.
//   ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//
//   auto specs = make_transform_specs_pt(pt, transform);
//   c_do_transform_pt(pt, specs);
// }
//
//
// // [[Rcpp::export]]
// void ParamTable_do_transform_premap(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
//   ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   // 1) Build full transform specs for all active parameters
//   auto full_specs = make_transform_specs_pt(pt, transform);
//
//   // 2) Filter specs to only premap trend parameters
//   auto premap_specs = filter_specs_by_param_set(pt, full_specs,
//                                                 te.plan.premap_trend_params);
//
//   // 3) Apply transforms in-place to those columns
//   c_do_transform_pt(pt, premap_specs);
// }
//
// // [[Rcpp::export]]
// void ParamTable_do_transform_pretransform(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
//   ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   // 1) Build full transform specs for all active parameters
//   auto full_specs = make_transform_specs_pt(pt, transform);
//
//   // 2) Filter specs to only premap trend parameters
//   auto pretransform_specs = filter_specs_by_param_set(pt, full_specs,
//                                                 te.plan.pretransform_trend_params);
//
//   // 3) Apply transforms in-place to those columns
//   c_do_transform_pt(pt, pretransform_specs);
// }
//
// // [[Rcpp::export]]
// void ParamTable_do_transform_posttransform(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
//   ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   // 1) Build full transform specs for all active parameters
//   auto full_specs = make_transform_specs_pt(pt, transform);
//
//   // 2) Filter specs to only premap trend parameters
//   auto posttransform_specs = filter_specs_by_param_set(pt, full_specs,
//                                                       te.plan.posttransform_trend_params);
//
//   // 3) Apply transforms in-place to those columns
//   c_do_transform_pt(pt, posttransform_specs);
// }
//
//
// // [[Rcpp::export]]
// void ParamTable_do_transform_postmap(SEXP param_table_ptr, SEXP trend_engine_ptr, Rcpp::List transform) {
//   ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   // 1) Build full transform specs for all active parameters
//   auto full_specs = make_transform_specs_pt(pt, transform);
//
//   // 2) Filter specs to only premap trend parameters
//   auto postmap_specs = complement_specs_for_premap(pt, full_specs,
//                                                 te.plan.premap_trend_params);
//
//
//   // 3) Apply transforms in-place to those columns
//   c_do_transform_pt(pt, postmap_specs);
// }
//
//
// // [[Rcpp::export]]
// void ParamTable_bind_trendops(SEXP param_table_ptr, SEXP trend_engine_ptr) {
//   ParamTable& pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   te.runtime.bind_all_ops_to_paramtable(pt);
// }
//
//
// // // [[Rcpp::export]]
// // Rcpp::List TrendEngine_run_premap_kernels_debug(SEXP param_table_ptr,
// //                                                 SEXP trend_engine_ptr) {
// //   ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
// //   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
// //
// //   const std::size_t n_ops = te.runtime.premap_ops.size();
// //   Rcpp::List out(n_ops);
// //   Rcpp::CharacterVector nms(n_ops);
// //
// //   for (std::size_t i = 0; i < n_ops; ++i) {
// //     TrendOpRuntime& op = te.runtime.premap_ops[i];
// //
// //     // ensure kernel is run
// //     te.runtime.run_kernels_for_op(op, pt);
// //
// //     const std::vector<double>& traj = op.kernel_ptr->get_output();
// //     Rcpp::NumericVector rtraj(traj.begin(), traj.end());
// //
// //     out[i] = rtraj;
// //     nms[i] = op.spec->target_param;
// //   }
// //
// //   out.attr("names") = nms;
// //   return out;
// // }
//
//
// // [[Rcpp::export]]
// void TrendEngine_apply_premap_bases(SEXP param_table_ptr,
//                                     SEXP trend_engine_ptr) {
//   ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   auto& runtime = te.runtime;
//   const std::size_t n_ops = runtime.premap_ops.size();
//
//   for (std::size_t i = 0; i < n_ops; ++i) {
//     TrendOpRuntime& op = runtime.premap_ops[i];
//     runtime.apply_base_for_op(op, pt);
//   }
// }
//
// // [[Rcpp::export]]
// void TrendEngine_apply_pretransform_bases(SEXP param_table_ptr,
//                                     SEXP trend_engine_ptr) {
//   ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   auto& runtime = te.runtime;
//   const std::size_t n_ops = runtime.pretransform_ops.size();
//
//   Rcpp::Rcout << "has_pretransform = " << runtime.has_pretransform()
//               << ", n_pretransform_ops = " << runtime.pretransform_ops.size() << "\n";
//
//   for (std::size_t i = 0; i < n_ops; ++i) {
//     TrendOpRuntime& op = runtime.pretransform_ops[i];
//     const TrendOpSpec& spec = *op.spec;
//     Rcpp::Rcout << "Pretransform op " << i
//                 << " target_param=" << spec.target_param
//                 << " base_type=" << spec.base_type
//                 << "\n";
//
//     runtime.apply_base_for_op(op, pt);
//   }
// }
//
// // [[Rcpp::export]]
// void TrendEngine_apply_posttransform_bases(SEXP param_table_ptr,
//                                           SEXP trend_engine_ptr) {
//   ParamTable&  pt = *Rcpp::XPtr<ParamTable>(param_table_ptr);
//   TrendEngine& te = *Rcpp::XPtr<TrendEngine>(trend_engine_ptr);
//
//   auto& runtime = te.runtime;
//   const std::size_t n_ops = runtime.posttransform_ops.size();
//
//   for (std::size_t i = 0; i < n_ops; ++i) {
//     TrendOpRuntime& op = runtime.posttransform_ops[i];
//     runtime.apply_base_for_op(op, pt);
//   }
// }
//
//
