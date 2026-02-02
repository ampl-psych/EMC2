#include <Rcpp>
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

inline std::vector<TransformSpec> make_transform_specs_from_full(
    NumericMatrix pars,
    CharacterVector full_names,
    const std::vector<TransformSpec>& full_specs)
{
  // Create a quick lookup name -> index in full_names/specs
  std::unordered_map<std::string,int> name_to_idx;
  for (int i = 0; i < full_names.size(); ++i) {
    name_to_idx[Rcpp::as<std::string>(full_names[i])] = i;
  }

  int ncol = pars.ncol();
  std::vector<TransformSpec> specs(ncol);
  CharacterVector cparnames = colnames(pars);
  for (int j = 0; j < ncol; j++) {
    std::string colname = Rcpp::as<std::string>(cparnames[j]);
    auto it = name_to_idx.find(colname);
    TransformSpec sp;
    sp.col_idx = j;
    if (it != name_to_idx.end()) {
      const TransformSpec& base = full_specs[it->second];
      sp.code  = base.code;
      sp.lower = base.lower;
      sp.upper = base.upper;
    } else {
      sp.code  = IDENTITY;
      sp.lower = 0.0;
      sp.upper = 1.0;
    }
    specs[j] = sp;
  }
  return specs;
}

std::vector<TransformSpec> make_transform_specs_for_paramtable(
    const ParamTable& pt,
    const Rcpp::List& transform)
{
  using namespace Rcpp;
  using std::string;

  // gather 'func', 'lower', 'upper' from transform list
  CharacterVector func_charvec = transform["func"];
  NumericVector lower_numvec   = transform["lower"];
  NumericVector upper_numvec   = transform["upper"];

  // Build name -> code
  std::unordered_map<string, TransformCode> codeMap;
  {
    CharacterVector fnames = func_charvec.names();
    for (int i = 0; i < func_charvec.size(); ++i) {
      string name = as<string>(fnames[i]);
      string f    = as<string>(func_charvec[i]);
      TransformCode code = IDENTITY;
      if (f == "exp")   code = EXP;
      else if (f == "pnorm") code = PNORM;
      codeMap.emplace(name, code);
    }
  }

  // Build name -> (lower, upper)
  std::unordered_map<string, std::pair<double,double>> boundMap;
  {
    CharacterVector ln = lower_numvec.names();
    for (int i = 0; i < lower_numvec.size(); ++i) {
      boundMap[ as<string>(ln[i]) ].first = lower_numvec[i];
    }
    CharacterVector un = upper_numvec.names();
    for (int i = 0; i < upper_numvec.size(); ++i) {
      boundMap[ as<string>(un[i]) ].second = upper_numvec[i];
    }
  }

  // Build specs only for currently active columns
  const int n_active = pt.n_params();
  std::vector<TransformSpec> specs;
  specs.reserve(n_active);

  for (int k = 0; k < n_active; ++k) {
    int base_idx = pt.active_cols[k];
    string colname = Rcpp::as<string>(pt.base_names[base_idx]);

    TransformSpec sp;
    sp.col_idx = base_idx;  // IMPORTANT: index into pt.base

    // code
    auto itc = codeMap.find(colname);
    sp.code = (itc != codeMap.end()) ? itc->second : IDENTITY;

    // bounds
    auto itb = boundMap.find(colname);
    if (itb != boundMap.end()) {
      sp.lower = itb->second.first;
      sp.upper = itb->second.second;
    } else {
      sp.lower = 0.0;  // your defaults
      sp.upper = 1.0;
    }

    specs.push_back(sp);
  }

  return specs;
}
std::vector<TransformSpec> make_transform_specs_for_paramtable_from_full(
    const ParamTable& pt,
    const Rcpp::CharacterVector& full_names,
    const std::vector<TransformSpec>& full_specs)
{
  using namespace Rcpp;
  using std::string;

  // Create name -> index into full_specs
  std::unordered_map<string,int> name_to_idx;
  name_to_idx.reserve(full_names.size());
  for (int i = 0; i < full_names.size(); ++i) {
    name_to_idx[ as<string>(full_names[i]) ] = i;
  }

  const int n_active = pt.n_params();
  std::vector<TransformSpec> specs;
  specs.reserve(n_active);

  for (int k = 0; k < n_active; ++k) {
    int base_idx = pt.active_cols[k];
    string colname = as<string>(pt.base_names[base_idx]);

    TransformSpec sp;
    sp.col_idx = base_idx;  // again: index into pt.base

    auto it = name_to_idx.find(colname);
    if (it != name_to_idx.end()) {
      const TransformSpec& base_sp = full_specs[it->second];
      sp.code  = base_sp.code;
      sp.lower = base_sp.lower;
      sp.upper = base_sp.upper;
    } else {
      sp.code  = IDENTITY;
      sp.lower = 0.0;
      sp.upper = 1.0;
    }

    specs.push_back(sp);
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
