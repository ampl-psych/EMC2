#include <Rcpp.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"
#include "trend.h"
#include "model_SS_EXG.h"
#include "race_integrate.h"
#include "model_SS_RDEX.h"
using namespace Rcpp;

LogicalVector c_do_bound(NumericMatrix pars,
                              const std::vector<BoundSpec>& specs)
{
  int nrows = pars.nrow();
  LogicalVector result(nrows, true);

  // For each parameter that has bounds
  for (size_t j = 0; j < specs.size(); j++) {
    const BoundSpec& bs = specs[j];
    int col_idx   = bs.col_idx;
    double min_v  = bs.min_val;
    double max_v  = bs.max_val;
    bool has_exc  = bs.has_exception;
    double exc_val= bs.exception_val;

    // Check each row
    for (int i = 0; i < nrows; i++) {
      double val = pars(i, col_idx);
      bool ok = (val > min_v && val < max_v);
      if (!ok && has_exc) {
        // If out of range, see if exception matches
        ok = (val == exc_val);
      }
      // Merge with existing result (like result = result & ok_col)
      if (result[i] && !ok) {
        result[i] = false;
      }
    }
  }
  return result;
}

NumericVector c_do_pre_transform(NumericVector p_vector,
                                 const std::vector<PreTransformSpec>& specs)
{
  for (size_t i = 0; i < specs.size(); i++) {
    const PreTransformSpec& s = specs[i];
    double val = p_vector[s.index];

    switch (s.code) {
    case PTF_EXP: {
      // lower + exp(real)
      p_vector[s.index] = s.lower + std::exp(val);
      break;
    }
    case PTF_PNORM: {
      double range = s.upper - s.lower;
      // lower + range * Φ(real)
      p_vector[s.index] = s.lower +
        range * R::pnorm(val, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
      break;
    }
    default:
      // no transform
      break;
    }
  }
  return p_vector;
}

NumericMatrix c_do_transform(NumericMatrix pars,
                             const std::vector<TransformSpec>& specs)
{
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
      // lower + exp(real)
      pars(i, col_idx) = lw + std::exp(pars(i, col_idx));
    }
      break;
    }
    case PNORM: {
      double range = up - lw;
      for (int i = 0; i < nrow; i++) {
        // lower + range * Φ(real)
        pars(i, col_idx) = lw +
          range * R::pnorm(pars(i, col_idx), 0.0, 1.0,
                           /*lower_tail=*/1, /*log_p=*/0);
      }
      break;
    }
    case IDENTITY:
    default:
      // do nothing
      break;
    }
  }
  return pars;
}


NumericMatrix c_map_p(NumericVector p_vector,
                      CharacterVector p_types,
                      List designs,
                      int n_trials,
                      DataFrame data,
                      List trend,
                      List transforms) {

  // Extract information about trends
  bool has_trend = (trend.length() > 0); // or another condition
  bool premap = false;
  bool pretransform = false;
  CharacterVector trend_names;
  // If trend has these flags
  if (has_trend) {
    premap = trend.attr("premap");
    pretransform = trend.attr("pretransform");
    trend_names = trend.names();
  }
  NumericVector p_mult_design;
  int n_params = p_types.size();
  NumericMatrix pars(n_trials, n_params);
  colnames(pars) = p_types;
  NumericMatrix trend_pars;
  // Identify trend parameters if any
  CharacterVector trend_pnames;
  LogicalVector trend_index(n_params, FALSE);
  if (has_trend && (premap || pretransform)) {
    // first loop over trends to get all trend pnames
    // But only for trends that are premap or pretransform
    for(unsigned int q = 0; q < trend.length(); q++){
      List cur_trend = trend[q];
      trend_pnames = c_add_charvectors(trend_pnames, as<CharacterVector>(cur_trend["trend_pnames"]));
      // Takes care of shared parameters
      trend_pnames = unique(trend_pnames);
    }
    // index which p_types are trends
    LogicalVector trend_index = contains_multiple(p_types,trend_pnames);
    for(unsigned int j = 0; j < trend_index.length(); j ++){
      // If we are a trend parameter:
      if(trend_index[j] == TRUE){
        NumericMatrix cur_design_trend = designs[j];
        CharacterVector cur_names_trend = colnames(cur_design_trend);
        // Take the current design and loop over columns
        // Multiply by design matrix
        for(int k = 0; k < cur_design_trend.ncol(); k ++){
          String cur_name_trend(cur_names_trend[k]);
          p_mult_design =  p_vector[cur_name_trend] * cur_design_trend(_, k);
          p_mult_design[is_nan(p_mult_design)] = 0;
          pars(_, j) = pars(_, j) + p_mult_design;
        }
      }
    }
    trend_pars = submat_rcpp_col_by_names(pars, trend_pnames);
    std::vector<TransformSpec> t_specs = make_transform_specs(trend_pars, transforms);
    trend_pars = c_do_transform(trend_pars, t_specs);
    trend_index = contains_multiple(p_types, trend_pnames);
  }
  for(int i = 0, t = 0; i < n_params; i++){
    if(trend_index[i] == FALSE){
      NumericMatrix cur_design = designs[i];
      CharacterVector cur_names = colnames(cur_design);

      for(int j = 0; j < cur_design.ncol(); j ++){
        String cur_name(cur_names[j]);
        NumericVector p_mult_design(n_trials, p_vector[cur_name]);
        // at this point we're multiplying by specific parameters (e.g. v_lMd)
        // So first apply trend to this parameter, then multiply by design matrix;
        if(has_trend && premap){
          // Check if trend is on current parameter
          LogicalVector cur_has_trend = contains(trend_names, cur_name);
          // This is a bit tricky and arguable.
          // Here we first fill a p_mult_design vector, then apply a trend then multiply with design matrix
          // Arguably you could also multiply parameter with design matrix and then apply trend
          // But that results in weird effects that if a parameter is set at 0, it could no longer be 0 post-trend
          for(unsigned int w = 0; w < cur_has_trend.length(); w ++){
            if(cur_has_trend[w] == TRUE){ // if so apply trend
              List cur_trend = trend[cur_name];
              CharacterVector cur_trend_pnames = cur_trend["trend_pnames"];
              p_mult_design = run_trend_rcpp(data, cur_trend, p_mult_design,
                                             submat_rcpp_col_by_names(trend_pars,cur_trend_pnames));

            }
          }
        }
        p_mult_design = p_mult_design * cur_design(_, j);
        p_mult_design[is_nan(p_mult_design)] = 0;
        pars(_, i) = pars(_, i) + p_mult_design;
      };
    } else if(pretransform){
      // These trends aren't applied here, but rather after mapping,
      // But they are transformed here already, so input them here.
      pars(_, i) = trend_pars(_, t);
      t++;
    }
  };
  if(has_trend && premap){
    pars = submat_rcpp_col(pars, !contains_multiple(p_types, trend_pnames));
  }
  return(pars);
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, List transforms, const std::vector<PreTransformSpec>& p_specs,
                              CharacterVector p_types, List designs, int n_trials, DataFrame data, List trend){
  bool has_trend = (trend.length() > 0);
  bool pretransform = false;
  bool posttransform = false;
  // If trend has these flags
  if (has_trend) {
    pretransform = trend.attr("pretransform");
    posttransform = trend.attr("posttransform");
  }
  NumericVector p_vector_updtd(clone(p_vector));
  CharacterVector par_names = p_vector_updtd.names();
  p_vector_updtd = c_do_pre_transform(p_vector_updtd, p_specs);
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, data, trend, transforms);
  // // Check if pretransform trend applies
  if(pretransform){ // automatically only applies if trend
    pars = prep_trend(data, trend, pars);
  }
  std::vector<TransformSpec> t_specs = make_transform_specs(pars, transforms);
  pars = c_do_transform(pars, t_specs);
  // Check if posttransform trend applies
  if(posttransform){ // automatically only applies if trend
    pars = prep_trend(data, trend, pars);
  }
  // ok is calculated afterwards and Ttransform applied in the function
  return(pars);
}

double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                            const int n_trials, IntegerVector expand,
                            double min_ll, LogicalVector is_ok){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  NumericVector lls(n_trials);
  NumericVector lls_exp(n_out);
  lls = d_DDM_Wien(rts, R, pars, is_ok);
  lls_exp = c_expand(lls, expand); // decompress
  // lls_exp = lls;
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return(sum(lls_exp));
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
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll), is_ok)); //first for compressed
  lds[winner] = win;
  if(n_acc > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll), is_ok)); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }
  lds[is_na(lds)] = min_ll;

  if(n_acc > 1){
    // LogicalVector winner_exp = c_bool_expand(winner, expand);
    NumericVector ll_out = lds[winner];
    NumericVector lds_los = lds[!winner];
    if(n_acc == 2){
      ll_out = ll_out + lds_los;
    } else{
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
      }
    }

    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    ll_out = c_expand(ll_out, expand); // decompress
    return(sum(ll_out));
  } else{
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
    lds_exp = c_expand(lds, expand); // decompress
    return(sum(lds_exp));
  }
}

double c_log_likelihood_ss(
    String type,
    NumericMatrix pars,
    DataFrame data,
    const int n_trials,
    IntegerVector expand,
    double min_ll,
    LogicalVector is_ok
) {
  // initialise local variables
  const int n_out = expand.length();
  if (is_true(all(!is_ok))) {
    NumericVector lls_expanded(n_out, min_ll);
    return(sum(lls_expanded));
  }
  NumericVector lls(n_trials);
  NumericVector lls_expanded(n_out);
  // extract data
  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  NumericVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  if (data.containsElementNamed("lI")) {
    IntegerVector lI = data["lI"];
    if (unique(lI).size() > 1) {
      stop("Column 'lI' contains multiple unique values; not yet supported.");
    }
  }
  // compute log likelihoods
  if (type == "SS_EXG") {
    lls = ss_exg_lpdf(RT, R, SSD, lR, winner, pars, is_ok, min_ll);
  } else if (type == "SS_TEXG") {
    lls = ss_texg_lpdf(RT, R, SSD, lR, winner, pars, is_ok, min_ll);
  } else {
    lls = ss_rdex_lpdf(RT, R, SSD, lR, winner, pars, is_ok, min_ll);
  }
  // decompress
  lls_expanded = c_expand(lls, expand);
  // protect against numerical issues
  lls_expanded = check_ll(lls_expanded, min_ll);

  // TODO REMOVE IF WE'RE SATISFIED WITH TESTING
//  Environment global = Environment::global_env();
//  global["llscpp"] = lls;


  // return summed log-likelihood
  return(sum(lls_expanded));
}



//Jeroen
double c_log_likelihood_race_missing(NumericMatrix pars, DataFrame data,
  NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
  NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
  const int n_trials, LogicalVector winner, IntegerVector expand,
  double min_ll, LogicalVector is_ok){


//  Rcpp::Rcout << "Entering missing" << std::endl;

  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  CharacterVector R = data["R"];
  NumericVector pCont = pars(_ , pars.ncol() - 1);
  NumericVector lds_exp(n_out);
  int n_acc = unique(R).length();
  if (any(is_na(R))) {
    n_acc -= 1;
  }
  if(sum(contains(data.names(), "NACC")) == 1){
    NumericVector lR = data["lR"];
    NumericVector NACC = data["NACC"];
    for(int x = 0; x < pars.nrow(); x++){
      if(lR[x] > NACC[x]){
        pars(x,0) = NA_REAL;
      }
    }
  }
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll), is_ok)); //first for compressed
  lds[winner] = win;
  if(n_acc > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll), is_ok)); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }

  lds[is_na(lds) | (winner & is_infinite(rts))] = min_ll;
  lds[(!winner) & (is_infinite(rts) | is_na(rts))] = 0;

//  Rcpp::Rcout << "lds 1 " << lds.length() << " " << lds << std::endl;

  // Calculate truncation?
  double LT = 0;
  double UT = R_PosInf;

  Nullable<double> LT_ = data.attr("LT");
  Nullable<double> UT_ = data.attr("UT");
  bool dotrunc = LT_.isNotNull() | UT_.isNotNull();

  if (LT_.isNotNull()) {
    LT = as<double>(LT_);
  }
  if (UT_.isNotNull()) {
    UT = as<double>(UT_);
  }

  // Calculate censoring
  double LC;
  double UC;
  Nullable<double> LC_ = data.attr("LC");
  Nullable<double> UC_ = data.attr("UC");

  if (LC_.isNotNull()) {
    LC = as<double>(LC_);
  }
  if (UC_.isNotNull()) {
    UC = as<double>(UC_);
  }

  // Response known
  // Fast
  LogicalVector neginf = rts == R_NegInf; // also sets NA to FALSE
  LogicalVector nortfast = neginf & !is_na(R);

  if (is_true(any(nortfast))) {

//    Rcpp::Rcout << "Response known fast" << std::endl;

    NumericMatrix mparsfast(sum(nortfast),pars.ncol());
    for (int i = 0, j = 0; i < nortfast.length(); i++) {
      if (nortfast[i]) {
        mparsfast(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnerfastvec = winner[nortfast];
    LogicalMatrix winnerfast(n_acc, winnerfastvec.length() / n_acc, winnerfastvec.begin());

    LogicalVector tofixfast = (winner & nortfast);
    NumericVector ldstofixfast(sum(tofixfast));

    for (int i = 0; i < sum(tofixfast); i++) {
      NumericMatrix pifast(n_acc, pars.ncol());
      if (n_acc == 1) {
        pifast(0,_) = mparsfast(i,_);
      } else {
        for (int j = 0; j < n_acc; j++) {
          pifast(j,_) = mparsfast(i * n_acc + j,_);
        }
      }
      NumericVector tmp = f_integrate(pifast, winnerfast(_,i), dfun, pfun, min_ll, LT, LC, is_ok);
      ldstofixfast[i] = std::log(std::max(0.0, std::min(tmp[0], 1.0)));
    }
    lds[tofixfast] = ldstofixfast;
  }

  // Slow
  LogicalVector posinf = rts == R_PosInf; // also sets NA to FALSE
  LogicalVector nortslow = posinf & !is_na(R);

  if (is_true(any(nortslow))) {

//    Rcpp::Rcout << "Response known slow" << std::endl;

    NumericMatrix mparsslow(sum(nortslow),pars.ncol());
    for (int i = 0, j = 0; i < nortslow.length(); i++) {
      if (nortslow[i]) {
        mparsslow(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnerslowvec = winner[nortslow];
    LogicalMatrix winnerslow(n_acc, winnerslowvec.length() / n_acc, winnerslowvec.begin());

    LogicalVector tofixslow = (winner & nortslow);
    NumericVector ldstofixslow(sum(tofixslow));

    for (int i = 0; i < sum(tofixslow); i++) {
      NumericMatrix pislow(n_acc, pars.ncol());
      if (n_acc == 1) {
        pislow(0,_) = mparsslow(i,_);
      } else {
        for (int j = 0; j < n_acc; j++) {
          pislow(j,_) = mparsslow(i * n_acc + j,_);
        }
      }
      NumericVector tmp = f_integrate_slow(pislow, winnerslow(_,i), dfun, pfun, min_ll, UC, UT, is_ok);
      ldstofixslow[i] = std::log(std::max(0.0, std::min(tmp[0], 1.0)));
    }
    lds[tofixslow] = ldstofixslow;
  }

  // No direction
  LogicalVector nortno = is_na(rts) & !is_na(R);

  if (is_true(any(nortno))) {

//    Rcpp::Rcout << "Response known no direction" << std::endl;

    NumericMatrix mparsno(sum(nortno), pars.ncol());

    for (int i = 0, j = 0; i < nortno.length(); i++) {
      if (nortno[i]) {
        mparsno(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnernovec = winner[nortno];
    LogicalMatrix winnerno(n_acc, winnernovec.length() / n_acc, winnernovec.begin());

    LogicalVector tofixno = (winner & nortno);
    NumericVector ldstofixno(sum(tofixno));

    for (int i = 0; i < sum(tofixno); i++) {
      NumericMatrix pino(n_acc, pars.ncol());
      if (n_acc  == 1) {
        pino(0,_) = mparsno(i,_);
      } else {
        for (int j = 0; j < n_acc; j++) {
          pino(j,_) = mparsno(i * n_acc + j,_);
        }
      }
      NumericVector tmpslow = f_integrate(pino, winnerno(_,i), dfun, pfun, min_ll, UC, UT, is_ok);
      NumericVector tmpfast = f_integrate(pino, winnerno(_,i), dfun, pfun, min_ll, LT, LC, is_ok);
      double tmp = tmpslow[0] + tmpfast[0];

      ldstofixno[i] = std::log(std::max(0.0, std::min(tmp, 1.0)));
    }
    lds[tofixno] = ldstofixno;
  }

  // Response unknown
  // Fast
  LogicalVector nortfastu = (rts == R_NegInf) & is_na(R);
  LogicalVector tofixfast(winner.length());

  if (is_true(any(nortfastu))) {

//    Rcpp::Rcout << "Response unknown fast" << std::endl;

    NumericMatrix mpars(sum(nortfastu), pars.ncol());
    for (int i = 0, j = 0; i < nortfastu.length(); i++) {
      if (nortfastu[i]) {
        mpars(j,_) = pars(i,_);
        j++;
      }
    }

//    Rcpp::Rcout << "dim pars " << pars.nrow() << " " << pars.ncol() << std::endl;
//    Rcpp::Rcout << "dim mpars " << mpars.nrow() << " " << mpars.ncol() << std::endl;
//    Rcpp::Rcout << "nortfastu " << nortfastu.length() << " " << nortfastu << std::endl;

    LogicalVector winnerfastuvec = winner[nortfastu];
    LogicalMatrix winnerfastu(n_acc, winnerfastuvec.length() / n_acc, winnerfastuvec.begin());

    tofixfast = (winner & nortfastu);
    NumericVector ldstofixfast(sum(tofixfast));

    for (int i = 0; i < sum(tofixfast); i++) {

//      NumericMatrix pi(n_acc, pars.nrow());
      NumericMatrix pi(n_acc, pars.ncol());

//      Rcpp::Rcout << n_acc << std::endl;
//      Rcpp::Rcout << pars.nrow() << std::endl;

//      Rcpp::Rcout << mpars.nrow() << std::endl;
//      Rcpp::Rcout << mpars.ncol() << std::endl;
//      Rcpp::Rcout << mpars << std::endl;


      for (int j = 0; j < n_acc; j++) {

//        Rcpp::Rcout << j << std::endl;

        pi(j,_) = mpars(i * n_acc + j,_);
      }

      LogicalVector idx(n_acc, 0);
      idx[0] = 1;

      NumericVector pc = f_integrate(pi, idx, dfun, pfun, min_ll, LT, LC, is_ok);

//      Rcpp::Rcout << "pc " << pc.length() << " " << pc << std::endl;

      double p;
//      if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
      if (pc[2] > 2 || traits::is_nan<REALSXP>(pc[0])) {
        p = NA_REAL;
      } else{
        p = std::max(0.0 ,std::min(pc[0],1.0));
      }

//      Rcpp::Rcout << "p " <<  p << std::endl;

      double cf;
      if (p != 0 && !(LT==0 && UT==R_PosInf)) {
        cf = pr_pt(pi, idx, dfun, pfun, min_ll, LT, UT, is_ok);
      } else {
        cf = 1;
      }

//      Rcpp::Rcout << "cf " << cf << std::endl;

      if (!traits::is_na<REALSXP>(cf)) {
        p *= cf;
      }

      if (!traits::is_na<REALSXP>(p) && n_acc > 1) {

        for (int j = 1; j < n_acc; j++) {

//          Rcpp::Rcout << "in loop " << std::endl;

          idx.fill(0);
          idx[j] = 1;
          pc = f_integrate(pi, idx, dfun, pfun, min_ll, LT, LC, is_ok);

//          Rcpp::Rcout << "pc " << pc.length() << " " << pc << std::endl;

          if (pc[2] > 2 || traits::is_nan<REALSXP>(pc[0])) {
            p = NA_REAL;
            break;
          }
          if (pc[0] != 0.0 && !(LT == 0 && UT == R_PosInf)) {
            cf = pr_pt(pi, idx, dfun, pfun, min_ll, LT, UT, is_ok);
          } else{
            cf = 1;
          }

//          Rcpp::Rcout << "cf " << cf << std::endl;

          if (!traits::is_na<REALSXP>(cf)) {
            p += pc[0] * cf;
          }
        }
      }
      double lp = std::log(p);
      if (!traits::is_na<REALSXP>(lp)) {
        ldstofixfast[i] = lp;
      } else{
        ldstofixfast[i] = R_NegInf;
      }
    }

//    Rcpp::Rcout << "tofixfast " << tofixfast.length() << " " << tofixfast << std::endl;
//    Rcpp::Rcout << "ldstofixfast " << ldstofixfast.length() << " " << ldstofixfast << std::endl;

   lds[tofixfast] = ldstofixfast;

//   Rcpp::Rcout << "lds 2 " << lds.length() << " " << lds << std::endl;

  }

  // Slow
  LogicalVector nortslowu = (rts == R_PosInf) & is_na(R);
  LogicalVector tofixslow(winner.length());

  if (is_true(any(nortslowu))) {

//    Rcpp::Rcout << "Response unknown slow" << std::endl;

    NumericMatrix mpars(sum(nortslowu), pars.ncol());

    for (int i = 0, j = 0; i < nortslowu.length(); i++) {
      if (nortslowu[i]) {
        mpars(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnerslowuvec = winner[nortslowu];
    LogicalMatrix winnerslowu(n_acc, winnerslowuvec.length() / n_acc, winnerslowuvec.begin());

    tofixslow = (winner & nortslowu);
    NumericVector ldstofixslow(sum(tofixslow));

    for (int i = 0; i < sum(tofixslow); i++) {
//      NumericMatrix pi(n_acc, pars.nrow());
      NumericMatrix pi(n_acc, pars.ncol());

      for (int j = 0; j < n_acc; j++) {
        pi(j,_) = mpars(i * n_acc + j,_);
      }

      LogicalVector idx(n_acc);
      idx[0] = 1;

      NumericVector pc = f_integrate(pi, idx, dfun, pfun, min_ll, UC, UT, is_ok);
      double p;
      if (pc[2] > 2 || traits::is_nan<REALSXP>(pc[0])) {
        p = NA_REAL;
      } else{
        p = std::max(0.0,std::min(pc[0],1.0));
      }

      double cf;
      if (p != 0 && !(LT==0 && UT==R_PosInf)) {
        cf = pr_pt(pi, idx, dfun, pfun, min_ll, LT, UT, is_ok);
      } else {
        cf = 1;
      }

      if (!traits::is_na<REALSXP>(cf)) {
        p *= cf;
      }

      if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
        for (int j = 1; j < n_acc; j++) {
          idx.fill(0);
          idx[j] = 1;
          pc = f_integrate(pi, idx, dfun, pfun, min_ll, UC, UT, is_ok);
          if (pc[2] > 2 || traits::is_nan<REALSXP>(pc[0])) {
            p = NA_REAL;
            break;
          }
          if (pc[0] != 0.0 && !(LT == 0 && UT == R_PosInf)) {
            cf = pr_pt(pi, idx, dfun, pfun, min_ll, LT, UT, is_ok);
          } else{
            cf = 1;
          }
          if (!traits::is_na<REALSXP>(cf)) {
            p += pc[0] * cf;
          }
        }
      }
      double lp = std::log(p);
      if (!traits::is_na<REALSXP>(lp)) {
        ldstofixslow[i] = lp;
      } else{
        ldstofixslow[i] = R_NegInf;
      }
    }
    lds[tofixslow] = ldstofixslow;
  }

  // No direction
  LogicalVector nortnou = is_na(rts) & is_na(R);
  nortnou = nortnou & (pCont == 0); // Otherwise non-identifiable

  if (is_true(any(nortnou))) {

//    Rcpp::Rcout << "Response unknown no direction" << std::endl;

    NumericMatrix mpars(sum(nortnou), pars.ncol());

    for (int i = 0, j = 0; i < nortnou.length(); i++) {
      if (nortnou[i]) {
        mpars(j,_) = pars(i,_);
        j++;
      }
    }

    LogicalVector winnernouvec = winner[nortnou];
    LogicalMatrix winnernou(n_acc, winnernouvec.length() / n_acc, winnernouvec.begin());

    LogicalVector tofix = (winner & nortnou);
    NumericVector ldstofix(sum(tofix));

    for (int i = 0; i < sum(tofix); i++) {
      NumericMatrix pi(n_acc, pars.ncol());

      for (int j = 0; j < n_acc; j++) {
        pi(j,_) = mpars(i * n_acc + j,_);
      }

      LogicalVector idx(n_acc);
      idx[0] = 1;

      double pc = pLU(pi, idx, dfun, pfun, min_ll, LT, LC, UC, UT, is_ok);
      double p;
      double cf;
      if (traits::is_na<REALSXP>(pc)) {
        p = NA_REAL;
      } else{
        if (pc != 0.0 && !(LT == 0 && UT == R_PosInf)) {
          cf = pr_pt(pi, idx, dfun, pfun, min_ll, LT, UT, is_ok);
        } else{
          cf = 1;
        }

        if (!traits::is_na<REALSXP>(cf)) {
          p = pc*cf;
        } else {
          p = NA_REAL;
        }

        if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
          for (int j = 1; j < n_acc; j++) {
            idx.fill(0);
            idx[j] = 1;
            pc = pLU(pi, idx, dfun, pfun, min_ll, LT, LC, UC, UT, is_ok);
            if (traits::is_na<REALSXP>(pc)) {
              p = NA_REAL;
              break;
            }
            if (pc != 0 && !(LT == 0 && UT == R_PosInf)) {
              cf = pr_pt(pi, idx, dfun, pfun, min_ll, LT, UT, is_ok);
            } else {
              cf = 1;
            }

            if (traits::is_na<REALSXP>(cf)) {
              p = NA_REAL;
              break;
            } else{
              p += pc * cf;
            }
          }
        }
      }
      double lp = std::log(p);
      if (!traits::is_na<REALSXP>(lp)) {
        ldstofix[i] = lp;
      } else{
        ldstofix[i] = R_NegInf;
      }
    }
    lds[tofix] = ldstofix;
  }

  // Truncation where not censored or censored and response known
  LogicalVector unique_nort = data.attr("unique_nort");
  NumericVector uniquewinlike = lds[unique_nort & winner];
  LogicalVector ok = is_finite(uniquewinlike);

//  Rcpp::Rcout << "unique_nort " << unique_nort.length() << " " << unique_nort << std::endl;
//  Rcpp::Rcout << "winner " << winner.length() << " " << winner << std::endl;
//  Rcpp::Rcout << "R " << R.length() << " " << R << std::endl;
//  Rcpp::Rcout << "ok " << ok.length() << " " << ok << std::endl;

  CharacterVector uniquewinresp = R[unique_nort & winner];
  LogicalVector alreadyfixed = is_na(uniquewinresp);
  ok = ok & !alreadyfixed;

//  Rcpp::Rcout << "uniquewinlike " << uniquewinlike.length() << " " << uniquewinlike << std::endl;
//  Rcpp::Rcout << "ok " << ok.length() << " " << ok << std::endl;
//  Rcpp::Rcout << "uniquewinresp " << uniquewinresp.length() << " " << uniquewinresp << std::endl;

  if (dotrunc & is_true(any(ok))) {


//    Rcpp::Rcout << "Truncation" << std::endl;

    IntegerVector expand_nort = data.attr("expand_nort");
    NumericMatrix tpars(sum(unique_nort),pars.ncol());
    for (int i=0, j=0; i < unique_nort.length(); i++) {
      if (unique_nort[i]) {
        tpars(j,_) = pars(i,_);
        j++;
      }
    }
    LogicalVector winnertruncvec = winner[unique_nort];
    LogicalMatrix winnertrunc(n_acc, winnertruncvec.length() / n_acc, winnertruncvec.begin());
    NumericVector cf = rep(NA_REAL, ok.length());

    for (int i = 0; i < ok.length(); i++) {
      if (ok[i]) {
        NumericMatrix pi(n_acc, pars.ncol());

        for (int j = 0; j < n_acc; j++) {
          pi(j,_) = tpars(i * n_acc + j , _ );
        }

        cf[i] = pr_pt(pi, winnertrunc(_,i), dfun, pfun, min_ll, LT, UT, is_ok);
      }
    }
    NumericVector cf_log = rep_each(log(cf), n_acc);
    NumericVector cf_exp = c_expand_jeroen(cf_log, expand_nort);
    LogicalVector fix = winner & !is_na(cf_exp) & is_finite(cf_exp);
    if (is_true(any(fix))) {
      lds[fix] = lds[fix] + cf_exp[fix];
    }
    LogicalVector badfix = winner & (is_na(cf_exp) | is_infinite(cf_exp));
    if (all(!is_na(tofixfast))) {
      badfix = badfix & !tofixfast;
    }
    if (all(!is_na(tofixslow))) {
      badfix = badfix & !tofixslow;
    }
    if (is_true(any(badfix))) {
      lds[badfix] = R_NegInf;
    }
  }

  // Non-process (contaminant) miss.
  LogicalVector isPCont = (pCont > 0);
  if (is_true(any(isPCont))) {

    // Rcpp::Rcout << "Contamination" << std::endl;

    NumericVector p = exp(lds[winner]);
    NumericVector pc = pCont[winner];
    CharacterVector Rwin = R[winner];
    LogicalVector isMiss = is_na(Rwin);
    for (int i = 0; i < p.length(); i++) {
      if (isMiss[i]) {
        p[i] = pc[i] + (1  - pc[i]) * p[i];
      } else {
        p[i] = (1 - pc[i]) * p[i];
      }
    }
    NumericVector ldswinner = log(p);
    lds[winner] = ldswinner;
  }

  // original code

// Rcpp::Rcout << "lds 3 = " << lds.length() << " " << lds << std::endl;

//  Rcpp::Rcout << "expand = " << expand.length() << " " << expand << std::endl;

//  lds_exp = c_expand_jeroen(lds, expand); // decompress

//  Rcpp::Rcout << "lds_exp 1 = " << lds_exp.length() << " " << lds_exp << std::endl;
//  Rcpp::Rcout << "n_acc  = " << n_acc << std::endl;

  if(n_acc > 1){

//    Rcpp::Rcout << "winner = " << winner.length() << " " << winner << std::endl;

    NumericVector ll_out = lds[winner];
    if(n_acc == 2){
      NumericVector lds_los = lds[!winner];
      ll_out = ll_out + lds_los;
    } else{
      NumericVector lds_los = lds_exp[!winner];
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
      }
    }

    // Rcpp::Rcout << "ll_out = " << ll_out.length() << " " << ll_out << std::endl;

    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;

//    Rcpp::Rcout << "ll_out = " << ll_out.length() << " " << ll_out << std::endl;

    lds_exp = c_expand_jeroen(ll_out, expand);

//    Rcpp::Rcout << "ll_exp final = " << lds_exp.length() << " " << lds_exp << std::endl;

    return(sum(lds_exp));
  } else{
    lds[is_na(lds)] = min_ll;
    lds[is_infinite(lds)] = min_ll;
    lds[lds < min_ll] = min_ll;
    lds_exp = c_expand_jeroen(lds, expand);
    return(sum(lds_exp));
  }


}


// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
            List designs, String type, List bounds, List transforms, List pretransforms,
            CharacterVector p_types, double min_ll, List trend){
  const int n_particles = p_matrix.nrow();
  const int n_trials = data.nrow();
  NumericVector lls(n_particles);
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  p_vector.names() = p_names;
  NumericMatrix pars(n_trials, p_types.length());
  LogicalVector is_ok(n_trials);

  // Once (outside the main loop over particles):
  NumericMatrix minmax = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<PreTransformSpec> p_specs;
  std::vector<BoundSpec> bound_specs;

  if (type == "DDM") {
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; i++) {
      p_vector = p_matrix(i, _);
      if (i == 0) {
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
      // Precompute specs
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if (type == "MRI" || type == "MRI_AR1") {
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for (int i = 0; i < n_particles; i++) {
      p_vector = p_matrix(i, _);
      if (i == 0) {
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      if (type == "MRI") {
        lls[i] = c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
      } else {
        lls[i] = c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
      }
    }
  } else if (type == "SS_TEXG" || type == "SS_EXG" || type == "SS_RDEX") {
    IntegerVector expand = data.attr("expand");
    NumericVector lR = data["lR"];
    const int n_lR = unique(lR).length();
    const int n_trials_ll = n_trials / n_lR;
    for (int i = 0; i < n_particles; i++) {
      p_vector = p_matrix(i, _);
      if (i == 0) {
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax, mm_names, pars, bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      lls[i] = c_log_likelihood_ss(type, pars, data, n_trials_ll, expand, min_ll, is_ok);
    }
  } else {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    // Love me some good old ugly but fast c++ pointers
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    if (type == "LBA" || type == "MLBA") {
      dfun = dlba_c;
      pfun = plba_c;
    } else if (type == "RDM" || type == "MRDM") {
      dfun = drdm_c;
      pfun = prdm_c;
    } else { // type == "LNR" || type == "MLNR"
      dfun = dlnr_c;
      pfun = plnr_c;
    }
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    for (int i = 0; i < n_particles; ++i) {
      p_vector = p_matrix(i, _);
      if (i == 0) {
        p_specs = make_pretransform_specs(p_vector, pretransforms);
      }
      pars = get_pars_matrix(p_vector, constants, transforms, p_specs, p_types, designs, n_trials, data, trend);
      if (i == 0) {
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      if (type == "MLBA" || type == "MRDM" || type == "MLNR" ) {
        lls[i] = c_log_likelihood_race_missing(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
      } else {
        lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
      }
    }
  }
  return(lls);
}

