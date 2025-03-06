#ifndef utility_h
#define utility_h

#include <Rcpp.h>
#include <unordered_set>
#include <vector>
#include <functional>
using namespace Rcpp;

LogicalVector contains(CharacterVector sv, std::string txt) {
  LogicalVector res(sv.size());
  for (int i = 0; i < sv.size(); i ++) {
    res[i] = (sv[i] == txt);
  }
  return res;
}

NumericVector vector_pow(NumericVector x1, NumericVector x2){
  NumericVector out(x1.length());
  for(unsigned int i = 0; i < out.length(); i ++){
    out[i] = pow(x1[i], x2[i]);
  }
  return(out);
}

NumericVector pnorm_multiple(NumericVector x){
  NumericVector out(x.size());
  for(int i = 0; i < x.size(); i++){
    out[i] = R::pnorm(x[i], 0, 1, TRUE, FALSE);
  }
  return out;
}

LogicalVector contains_multiple(CharacterVector sv, CharacterVector inputs) {
  LogicalVector res(sv.size());
  for (int i = 0; i < sv.size(); i ++) {
    int k = 0;
    for (int j = 0; j < inputs.size(); j++){
      if (sv[i] == inputs[j]){
        k++;
      }
    }
    res[i] = k > 0;
  }
  return res;
}

NumericMatrix submat_rcpp_col(NumericMatrix X, LogicalVector condition) {
  int n = X.nrow();
  int k = X.ncol();

  if (condition.size() != k) {
    stop("Length of logical vector must match number of columns of the matrix.");
  }

  int to_keep = sum(condition);

  // If all columns match, just return X to avoid unnecessary copying
  if (to_keep == k) {
    return X;
  }

  NumericMatrix out(n, to_keep);

  double* x_ptr = REAL(X);
  double* out_ptr = REAL(out);

  // We'll keep track of the next column in 'out' to fill
  int out_col_index = 0;

  // Extract the original column names
  CharacterVector orig_colnames = colnames(X);
  CharacterVector new_colnames(to_keep);

  for (int col = 0; col < k; col++) {
    if (condition[col]) {
      double* x_col_start = x_ptr + col * n;
      double* out_col_start = out_ptr + out_col_index * n;

      std::copy(x_col_start, x_col_start + n, out_col_start);

      // Assign the column name
      new_colnames[out_col_index] = orig_colnames[col];

      out_col_index++;
    }
  }

  // Set new column names on output
  colnames(out) = new_colnames;

  return out;
}



NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition) {
  int n = X.nrow(), k = X.ncol();
  int to_keep = sum(condition);
  // If all rows match, just return X (this avoids copying)
  if (to_keep == n) {
    return X;
  }
  NumericMatrix out(to_keep, k);
  for (int i = 0, j = 0; i < n; i++) {
    if (condition[i]) {
      out(j, _) = X(i, _);
      j++;
    }
  }
  colnames(out) = colnames(X);
  return out;
}


NumericVector c_expand(NumericVector x1, IntegerVector expand){
  const int n_out = expand.length();
  NumericVector out(n_out);
  int curr_idx;
  for(int i = 0; i < n_out; i++){
    curr_idx = expand[i] - 1; //expand created in 1-based R
    out[i] = x1[curr_idx];
  }
  return(out);
}

LogicalVector c_bool_expand(LogicalVector x1, IntegerVector expand){
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

// [[Rcpp::export]]
CharacterVector c_add_charvectors(CharacterVector x, CharacterVector y) {
  // Create a new vector of length = length(x) + length(y)
  CharacterVector z(x.size() + y.size());
  // Copy x into z
  std::copy(x.begin(), x.end(), z.begin());
  // Copy y into z after x
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  return z;
}

// Custom hash for a vector<double>
struct RowHash {
  std::size_t operator()(const std::vector<double> &v) const {
    std::size_t seed = 0;
    std::hash<double> hash_double;
    for (double d : v) {
      // A standard hash combination approach (based on boost::hash_combine)
      seed ^= hash_double(d) + 0x9e3779b97f4a7c16ULL + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// Custom equality for a vector<double>
struct RowEqual {
  bool operator()(const std::vector<double> &a, const std::vector<double> &b) const {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
      if (a[i] != b[i]) return false;
    }
    return true;
  }
};

Rcpp::LogicalVector duplicated_matrix(Rcpp::NumericMatrix x) {
  int n = x.nrow();
  int m = x.ncol();

  Rcpp::LogicalVector dup(n);
  dup.fill(false);

  // Pre-allocate storage to reduce reallocations
  std::vector<double> buffer(m);
  std::unordered_set<std::vector<double>, RowHash, RowEqual> seen;
  seen.reserve(n);

  for (int i = 0; i < n; i++) {
    // Copy row data into buffer
    for (int j = 0; j < m; j++) {
      buffer[j] = x(i, j);
    }

    // Check if we have seen this row before
    if (seen.find(buffer) != seen.end()) {
      dup[i] = true;
    } else {
      // Insert a copy of the current row
      seen.insert(std::vector<double>(buffer.begin(), buffer.end()));
    }
  }

  return dup;
}

IntegerVector cumsum_logical(LogicalVector x) {
  int n = x.size();
  IntegerVector out(n);
  int running_total = 0;

  for (int i = 0; i < n; i++) {
    // Add 1 if TRUE, else 0
    if (x[i]) {
      running_total += 1;
    }
    out[i] = running_total;
  }

  return out;
}

IntegerVector which_rcpp(LogicalVector x) {
  int n = x.size();
  int count = 0;
  // First pass: count how many TRUE
  for (int i = 0; i < n; i++) {
    if (x[i]) count++;
  }

  // Allocate the output vector
  IntegerVector out(count);

  // Second pass: fill the output with the indices of TRUE values
  for (int i = 0, j = 0; i < n; i++) {
    if (x[i]) {
      out[j] = i;
      j++;
    }
  }

  return out;
}


// For do_bounds
struct BoundSpec {
  int col_idx;         // which column in 'pars'
  double min_val;
  double max_val;
  bool has_exception;
  double exception_val;
};


std::vector<BoundSpec> make_bound_specs(NumericMatrix minmax,
                                        CharacterVector minmax_colnames,
                                        NumericMatrix pars,
                                        List bound)
{
  // 1) Build a map from param-name -> column index in 'pars'
  CharacterVector pcolnames = colnames(pars);
  std::unordered_map<std::string, int> colMap;
  for (int j = 0; j < pcolnames.size(); j++) {
    colMap[ Rcpp::as<std::string>(pcolnames[j]) ] = j;
  }

  // 2) Build a map from param-name -> exception value
  bool has_exception = bound.containsElementNamed("exception") && !Rf_isNull(bound["exception"]);
  std::unordered_map<std::string, double> exceptionMap;
  if (has_exception) {
    NumericVector except_vec = bound["exception"];
    CharacterVector except_names = except_vec.names();
    for (int i = 0; i < (int)except_vec.size(); i++) {
      exceptionMap[ Rcpp::as<std::string>(except_names[i])] = except_vec[i];
    }
  }

  // 3) Create BoundSpec for each column in minmax
  int ncols = minmax_colnames.size();
  std::vector<BoundSpec> specs(ncols);
  for (int j = 0; j < ncols; j++) {
    std::string var_name = Rcpp::as<std::string>(minmax_colnames[j]);

    // Fill the struct
    BoundSpec s;
    s.col_idx     = colMap[var_name];
    s.min_val     = minmax(0, j);
    s.max_val     = minmax(1, j);

    auto it = exceptionMap.find(var_name);
    if (it != exceptionMap.end()) {
      s.has_exception = true;
      s.exception_val = it->second;
    } else {
      s.has_exception = false;
      s.exception_val = NA_REAL;  // or 0
    }
    specs[j] = s;
  }
  return specs;
}

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

std::vector<TransformSpec> make_transform_specs(NumericMatrix pars, List transform)
{
  // gather 'func', 'lower', 'upper'
  CharacterVector func_charvec = transform["func"];
  NumericVector lower_numvec   = transform["lower"];
  NumericVector upper_numvec   = transform["upper"];

  // Build a map param_name -> code
  std::unordered_map<std::string,TransformCode> codeMap;
  {
    // e.g. "param_name" -> "exp" or "pnorm" in func_charvec
    CharacterVector fnames = func_charvec.names();
    for (int i = 0; i < func_charvec.size(); i++) {
      std::string name = Rcpp::as<std::string>(fnames[i]);
      std::string f    = Rcpp::as<std::string>(func_charvec[i]);
      if (f == "exp") {
        codeMap[name] = EXP;
      } else if (f == "pnorm") {
        codeMap[name] = PNORM;
      } else {
        codeMap[name] = IDENTITY;
      }
    }
  }

  // Build param_name -> (lower, upper)
  std::unordered_map<std::string,std::pair<double,double>> boundMap;
  {
    CharacterVector ln = lower_numvec.names();
    for (int i = 0; i < lower_numvec.size(); i++) {
      std::string nm = Rcpp::as<std::string>(ln[i]);
      boundMap[nm].first = lower_numvec[i];
    }
    CharacterVector un = upper_numvec.names();
    for (int i = 0; i < upper_numvec.size(); i++) {
      std::string nm = Rcpp::as<std::string>(un[i]);
      boundMap[nm].second = upper_numvec[i];
    }
  }

  // Now fill specs for each col in pars
  int ncol = pars.ncol();
  std::vector<TransformSpec> specs(ncol);

  CharacterVector cparnames = colnames(pars);
  for (int j = 0; j < ncol; j++) {
    std::string colname = Rcpp::as<std::string>(cparnames[j]);
    TransformSpec sp;
    sp.col_idx = j;
    sp.code    = codeMap[colname];
    auto it    = boundMap.find(colname);
    if (it != boundMap.end()) {
      sp.lower = it->second.first;
      sp.upper = it->second.second;
    } else {
      sp.lower = 0.0;  // or default
      sp.upper = 1.0;  // or default
    }
    specs[j] = sp;
  }

  return specs;
}


// For pretransform
enum PreTFCode { PTF_EXP = 1, PTF_PNORM = 2, PTF_NONE = 0 };

struct PreTransformSpec {
  int index;       // index in p_vector
  PreTFCode code;
  double lower;
  double upper;
  // Possibly store the original name if needed
};

std::vector<PreTransformSpec> make_pretransform_specs(NumericVector p_vector, List transform)
{
  // e.g. transform["func"], transform["lower"], transform["upper"]
  CharacterVector func   = transform["func"];
  NumericVector lowervec = transform["lower"];
  NumericVector uppervec = transform["upper"];

  // Build a map param_name -> code
  std::unordered_map<std::string,PreTFCode> codeMap;
  CharacterVector fnames = func.names();
  for (int i = 0; i < func.size(); i++) {
    std::string name = Rcpp::as<std::string>(fnames[i]);
    std::string f    = Rcpp::as<std::string>(func[i]);
    if (f == "exp") {
      codeMap[name] = PTF_EXP;
    } else if (f == "pnorm") {
      codeMap[name] = PTF_PNORM;
    } else {
      codeMap[name] = PTF_NONE;
    }
  }

  // Build a map param_name -> (lower, upper)
  std::unordered_map<std::string, std::pair<double,double>> boundMap;
  {
    CharacterVector ln = lowervec.names();
    for (int i = 0; i < lowervec.size(); i++) {
      boundMap[ Rcpp::as<std::string>(ln[i]) ].first = lowervec[i];
    }
    CharacterVector un = uppervec.names();
    for (int i = 0; i < uppervec.size(); i++) {
      boundMap[ Rcpp::as<std::string>(un[i]) ].second = uppervec[i];
    }
  }

  // Now create PreTransformSpec for each element in p_vector
  CharacterVector p_names = p_vector.names();
  int n = p_vector.size();
  std::vector<PreTransformSpec> specs(n);
  for (int i = 0; i < n; i++) {
    std::string pname = Rcpp::as<std::string>(p_names[i]);
    PreTransformSpec s;
    s.index = i;
    s.code = codeMap[pname];
    auto it = boundMap.find(pname);
    if (it != boundMap.end()) {
      s.lower = it->second.first;
      s.upper = it->second.second;
    } else {
      s.lower = 0.0;
      s.upper = 1.0;
    }
    specs[i] = s;
  }
  return specs;
}



#endif


