#include <Rcpp.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "trend.h"
using namespace Rcpp;

LogicalVector c_do_bound(NumericMatrix pars, List bound) {
  // Extract 'minmax' matrix and its column names
  NumericMatrix minmax = bound["minmax"];
  CharacterVector minmax_colnames = colnames(minmax);

  // Map parameter names to column indices in 'pars'
  CharacterVector pars_colnames = colnames(pars);
  std::unordered_map<std::string, int> col_indices;
  for (int j = 0; j < pars_colnames.size(); ++j) {
    col_indices[as<std::string>(pars_colnames[j])] = j;
  }

  int nrows = pars.nrow();
  int ncols = minmax_colnames.size();

  // Initialize the result vector to 'true'
  LogicalVector result(nrows, true);

  // Extract exception vector and its names if 'exception' is not NULL
  bool has_exception = bound.containsElementNamed("exception") && !Rf_isNull(bound["exception"]);
  NumericVector exception_vec;
  CharacterVector exception_names;
  if (has_exception) {
    exception_vec = bound["exception"];
    exception_names = exception_vec.names();
  }

  // Loop over each parameter in 'minmax'
  for (int j = 0; j < ncols; ++j) {
    std::string var_name = as<std::string>(minmax_colnames[j]);
    int col_idx = col_indices[var_name];
    NumericVector pars_col = pars(_, col_idx);

    double min_val = minmax(0, j);
    double max_val = minmax(1, j);

    // Perform bounds checking
    LogicalVector ok_col = (pars_col > min_val) & (pars_col < max_val);

    // Apply exception if any
    if (has_exception && exception_vec.containsElementNamed(var_name.c_str())) {
      double exception_value = exception_vec[var_name];
      ok_col = ok_col | (pars_col == exception_value);
    }

    // Update the result vector
    result = result & ok_col;
  }

  return result;
}

NumericVector c_do_pre_transform(NumericVector p_vector,
                                 List transform) {
  // Extract transform components
  CharacterVector func = transform["func"];
  NumericVector lower = transform["lower"];
  NumericVector upper = transform["upper"];
  // Get the names of p_vector and transform components
  CharacterVector p_names = p_vector.names();
  CharacterVector func_names = func.names();
  CharacterVector lower_names = lower.names();
  CharacterVector upper_names = upper.names();

  int n = p_vector.size();
  // Match p_names to func_names to align indices
  // match() returns 1-based indices
  Function match("match");
  IntegerVector idx = match(p_names, func_names);

  if (is_true(any(is_na(idx)))) {
    stop("Some names in p_vector not found in transform$func.");
  }

  // Now, idx[i]-1 gives the position in func, lower, upper corresponding to p_vector[i]
  for (int i = 0; i < n; i++) {
    int pos = idx[i] - 1;
    std::string f = as<std::string>(func[pos]);
    if (f == "exp") {
      double val = p_vector[i];
      double lw = lower[pos];
      p_vector[i] = exp(val - lw);
    } else if (f == "pnorm") {
      double val = p_vector[i];
      double lw = lower[pos];
      double up = upper[pos];
      // Transform the value before pnorm
      double z = (val - lw) / (up - lw);
      p_vector[i] = R::pnorm(z, 0.0, 1.0, 1, 0);
    }
    // If other transformations exist, handle them similarly
  }
  return p_vector;
}

NumericMatrix c_do_transform(NumericMatrix pars, List transform) {
  // Get the column names of pars
  CharacterVector colnames_pars = colnames(pars);
  int ncol = pars.ncol();
  int nrow = pars.nrow();
  // Get 'func', 'lower', 'upper' from transform
  CharacterVector func_charvec = transform["func"];
  CharacterVector names_func = func_charvec.names();
  NumericVector lower_numvec = transform["lower"];
  CharacterVector names_lower = lower_numvec.names();
  NumericVector upper_numvec = transform["upper"];
  CharacterVector names_upper = upper_numvec.names();
  // Create maps from names to values
  std::map<String, String> func_map;
  for (int i = 0; i < func_charvec.size(); i++) {
    func_map[ names_func[i] ] = func_charvec[i];
  }
  std::map<String, double> lower_map;
  for (int i = 0; i < lower_numvec.size(); i++) {
    lower_map[ names_lower[i] ] = lower_numvec[i];
  }
  std::map<String, double> upper_map;
  for (int i = 0; i < upper_numvec.size(); i++) {
    upper_map[ names_upper[i] ] = upper_numvec[i];
  }
  // For each column, apply the appropriate transformation
  for (int j = 0; j < ncol; j++) {
    String colname = colnames_pars[j];
    String f = func_map[colname];
    double l = lower_map[colname];
    double u = upper_map[colname];

    if (f == "exp") {
      for (int i = 0; i < nrow; i++) {
        pars(i, j) = exp(pars(i,j) - l);
      }
    } else if (f == "pnorm") {
      double range = u - l;
      for (int i = 0; i < nrow; i++) {
        double z = (pars(i,j) - l) / range;
        pars(i,j) = R::pnorm(z, 0.0, 1.0, 1, 0);
      }
    }
    // For "identity" no transformation is applied
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
    // First deal with trend parameters
    for(unsigned int q = 0; q < trend.length(); q++){
      // Loop over trends
      List cur_trend = trend[q];
      LogicalVector cur_trend_idx = contains_multiple(p_types,as<CharacterVector>(cur_trend["trend_pnames"]));
      trend_pnames = c_add_charvectors(trend_pnames, as<CharacterVector>(cur_trend["trend_pnames"]));
      // Loop over p_types, pick out any that are trends
      for(unsigned int j = 0; j < cur_trend_idx.length(); j ++){
        if(cur_trend_idx[j] == TRUE){
          NumericMatrix cur_design_trend = designs[j];
          CharacterVector cur_names_trend = colnames(cur_design_trend);
          // Make an empty 1 column matrix
          // This is needed for c_do_transform
          // Take the current design and loop over columns
          for(int k = 0; k < cur_design_trend.ncol(); k ++){
            String cur_name_trend(cur_names_trend[k]);
            p_mult_design =  p_vector[cur_name_trend] * cur_design_trend(_, k);
            p_mult_design[is_nan(p_mult_design)] = 0;
            pars(_, j) = pars(_, j) + p_mult_design;
          }
        }
      }
    }
    trend_pars = c_do_transform(submat_rcpp_col(pars, contains_multiple(p_types, trend_pnames)), transforms);
    trend_index = contains_multiple(p_types, trend_pnames);
  }
  for(int i = 0, t = 0; i < n_params; i++){
    if(trend_index[i] == FALSE){
      NumericMatrix cur_design = designs[i];
      CharacterVector cur_names = colnames(cur_design);
      for(int j = 0; j < cur_design.ncol(); j ++){
        String cur_name(cur_names[j]);
        p_mult_design =  p_vector[cur_name] * cur_design(_, j);
        if(has_trend && premap){
          // Check if trend is on current parameter
          LogicalVector cur_has_trend = contains(trend_names, cur_name);
          for(unsigned int w = 0; w < cur_has_trend.length(); w ++){
            if(cur_has_trend[w] == TRUE){ // if so apply trend
              List cur_trend = trend[cur_name];
              CharacterVector cur_trend_pnames = cur_trend["trend_pnames"];
              p_mult_design = run_trend_rcpp(data, cur_trend, p_mult_design,
                                             submat_rcpp_col(trend_pars, contains_multiple(trend_pnames, cur_trend_pnames)));
            }
          }
        }
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
  if(has_trend){
    if(premap){
      pars = submat_rcpp_col(pars, !contains_multiple(p_types, trend_pnames));
    } else{

    }

  }
  return(pars);
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, List transforms, List pretransforms,
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
  p_vector_updtd = c_do_pre_transform(p_vector_updtd, pretransforms);
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, data, trend, transforms);
  // // Check if pretransform trend applies
  if(pretransform){ // automatically only applies if trend
    pars = prep_trend(data, trend, pars);
  }
  pars = c_do_transform(pars, transforms);
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
  NumericVector lds_exp(n_out);
  const int n_acc = unique(R).length();
  if(sum(contains(data.names(), "RACE")) == 1){
    NumericVector lR = data["lR"];
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
  lds_exp = c_expand(lds, expand); // decompress
  if(n_acc > 1){
    LogicalVector winner_exp = c_bool_expand(winner, expand);
    NumericVector ll_out = lds_exp[winner_exp];
    if(n_acc == 2){
      NumericVector lds_los = lds_exp[!winner_exp];
      ll_out = ll_out + lds_los;
    } else{
      NumericVector lds_los = lds_exp[!winner_exp];
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
      }
    }
    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    return(sum(ll_out));
  } else{
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
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
  IntegerVector expand = data.attr("expand");
  LogicalVector is_ok(n_trials);
  if(type == "DDM"){
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      pars = get_pars_matrix(p_vector, constants, transforms, pretransforms, p_types, designs, n_trials, data, trend);
      is_ok = c_do_bound(pars, bounds);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else{
    LogicalVector winner = data["winner"];
    // Love me some good old ugly but fast c++ pointers
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    // NumericMatrix (*Ttransform)(NumericMatrix);
    if(type == "LBA"){
      dfun = dlba_c;
      pfun = plba_c;
    } else if(type == "RDM"){
      dfun = drdm_c;
      pfun = prdm_c;
    } else{
      dfun = dlnr_c;
      pfun = plnr_c;
    }
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      pars = get_pars_matrix(p_vector, constants, transforms, pretransforms, p_types, designs, n_trials, data, trend);
      is_ok = c_do_bound(pars, bounds);
      lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
    }
  }
  return(lls);
}
