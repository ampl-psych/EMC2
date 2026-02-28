// transform_utils.h
#ifndef EMC2_TRANSFORM_UTILS_H
#define EMC2_TRANSFORM_UTILS_H

#include <Rcpp.h>
// #include <RcppArmadillo.h>
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

//

enum PreTFCode { PTF_EXP = 1, PTF_PNORM = 2, PTF_NONE = 0 };

struct PreTransformSpec {
  int index;      // index in p_vector
  PreTFCode code;
  double lower;
  double upper;
};

//

struct BoundSpec {
  int col_idx;
  double min_val;
  double max_val;
  bool has_exception;
  double exception_val;
};


// --- Declarations: matrix-based transforms ---

std::vector<TransformSpec> make_transform_specs(const Rcpp::NumericMatrix& pars,
                       const Rcpp::List& transform);

std::vector<TransformSpec> make_transform_specs_from_full(const Rcpp::NumericMatrix& pars,
                                 const Rcpp::CharacterVector& full_names,
                                 const std::vector<TransformSpec>& full_specs);

std::vector<PreTransformSpec> make_pretransform_specs(const Rcpp::NumericVector& p_vector,
                          const Rcpp::List& transform);

Rcpp::NumericVector c_do_pre_transform(Rcpp::NumericVector p_vector,
                   const std::vector<PreTransformSpec>& specs);

Rcpp::NumericMatrix c_do_transform(Rcpp::NumericMatrix pars,
               const std::vector<TransformSpec>& specs);

Rcpp::LogicalVector c_do_bound(Rcpp::NumericMatrix pars, const std::vector<BoundSpec>& specs);

Rcpp::LogicalVector c_do_bound_pt(const ParamTable& pt, const std::vector<BoundSpec>& specs);

// --- Declarations: ParamTable-based transforms ---

std::vector<TransformSpec>
  make_transform_specs_for_paramtable(const ParamTable& pt,
                                      const Rcpp::List& transform);

std::vector<TransformSpec>
  make_transform_specs_for_paramtable_from_full(
    const ParamTable& pt,
    const Rcpp::CharacterVector& full_names,
    const std::vector<TransformSpec>& full_specs);

void c_do_transform_pt(ParamTable& pt, const std::vector<TransformSpec>& specs);


// --- Declarations: Bound spec ---
std::vector<BoundSpec> make_bound_specs(Rcpp::NumericMatrix minmax,
                                        Rcpp::CharacterVector minmax_colnames,
                                        Rcpp::NumericMatrix pars,
                                        Rcpp::List bound);

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
