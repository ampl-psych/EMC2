// transform_utils.h
#ifndef EMC2_TRANSFORM_UTILS_H
#define EMC2_TRANSFORM_UTILS_H

#include <Rcpp.h>
#include "ParamTable.h"

// For transforms
enum TransformCode {
  IDENTITY = 0,
  EXP      = 1,
  PNORM    = 2
};

struct TransformSpec {
  int col_idx;        // which column in 'pars'
  TransformCode code; // e.g. EXP, PNORM, ...
  double lower;
  double upper;
};

struct BoundSpec {
  int col_idx;   // column in pt
  double min_val;
  double max_val;
  bool has_exception;
  double exception_val;
};

// --- Declarations: transforms ---
std::vector<TransformSpec> make_transform_specs_matrix(const Rcpp::NumericMatrix& pars,
                                                       const Rcpp::List& transform);

Rcpp::NumericMatrix c_do_transform_matrix(Rcpp::NumericMatrix pars,
                                          const std::vector<TransformSpec>& specs);

Rcpp::LogicalVector c_do_bound_pt(const ParamTable& pt, const std::vector<BoundSpec>& specs);

std::vector<TransformSpec> make_transform_specs_pt(const ParamTable& pt,
                                                   const Rcpp::List& transform);

void c_do_transform_pt(ParamTable& pt, const std::vector<TransformSpec>& specs);

std::vector<BoundSpec> make_bound_specs_pt(Rcpp::NumericMatrix minmax,
                                           Rcpp::CharacterVector minmax_colnames,
                                           const ParamTable& pt,
                                           Rcpp::List bound);

std::vector<TransformSpec> filter_specs_by_param_set(
    const ParamTable& pt,
    const std::vector<TransformSpec>& full_specs,
    const std::unordered_set<std::string>& allowed);

std::vector<TransformSpec> complement_specs_for_premap(
    const ParamTable& pt,
    const std::vector<TransformSpec>& full_specs,
    const std::unordered_set<std::string>& premap_set);


#endif
