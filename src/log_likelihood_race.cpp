#include <Rcpp.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
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
        pars(i,j) = exp(pars(i,j) - l);
      }
    } else if (f == "pnorm") {
      double range = u - l;
      for (int i = 0; i < nrow; i++) {
        double z = (pars(i,j) - l) / range;
        pars(i,j) = R::pnorm(z, 0.0, 1.0, 1, 0);
      }
    }
    // For "identity" or other functions, no transformation is applied
  }
  return pars;
}



NumericVector c_expand(NumericVector x1, NumericVector expand){
  const int n_out = expand.length();
  NumericVector out(n_out);
  int curr_idx;
  for(int i = 0; i < n_out; i++){
    curr_idx = expand[i] - 1; //expand created in 1-based R
    out[i] = x1[curr_idx];
  }
  return(out);
}

LogicalVector c_bool_expand(LogicalVector x1, NumericVector expand){
  const int n_out = expand.length();
  LogicalVector out(n_out);
  int curr_idx;
  for(int i = 0; i < n_out; i++){
    curr_idx = expand[i] - 1; //expand created in 1-based R
    out[i] = x1[curr_idx];
  }
  return(out);
}

NumericVector c_add_vectors(NumericVector x1, NumericVector x2){
  if(is_na(x2)[0] ){
    return(x1);
  }
  NumericVector output(x1.size() + x2.size());
  std::copy(x1.begin(), x1.end(), output.begin());
  std::copy(x2.begin(), x2.end(), output.begin() + x1.size());
  CharacterVector all_names(x1.size() + x2.size());
  CharacterVector x1_names = x1.names();
  CharacterVector x2_names = x2.names();
  std::copy(x1_names.begin(), x1_names.end(), all_names.begin());
  std::copy(x2_names.begin(), x2_names.end(), all_names.begin() + x1.size());
  output.names() = all_names;
  return output;
}

// LL generic functions
// [[Rcpp::export]]
NumericMatrix c_map_p(NumericVector p_vector, CharacterVector p_types, List designs, int n_trials){
  NumericMatrix pars(n_trials, p_types.length());
  for(int i = 0; i < p_types.length(); i++){
    NumericMatrix curr_design = designs[i];
    CharacterVector curr_names = colnames(curr_design);
    for(int j = 0; j < curr_design.ncol(); j ++){
      String curr_name(curr_names[j]);
      pars(_, i) = pars(_, i) + p_vector[curr_name] * curr_design(_, j);
    };
  };
  colnames(pars) = p_types;
  return(pars);
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, List transforms,
                       CharacterVector p_types, List designs, int n_trials){
  NumericVector p_vector_updtd(clone(p_vector));
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials);
  pars = c_do_transform(pars, transforms);
  return(pars);
}

double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                          const int n_trials, NumericVector expand,
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
                             const int n_trials, LogicalVector winner, NumericVector expand,
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
                      List designs, String type, List bounds, List transforms, CharacterVector p_types,
                      double min_ll, List group_idx){
  const int n_particles = p_matrix.nrow();
  const int n_trials = data.nrow();
  NumericVector lls(n_particles);
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  NumericMatrix pars(n_trials, p_types.length());
  p_vector.names() = p_names;
  NumericVector expand = data.attr("expand");
  LogicalVector is_ok(n_trials);
  if(type == "DDM"){
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      pars = get_pars_matrix(p_vector, constants, transforms, p_types, designs, n_trials);
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
      pars = get_pars_matrix(p_vector, constants, transforms, p_types, designs, n_trials);
      is_ok = c_do_bound(pars, bounds);
      lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
    }
  }
  return(lls);
}

