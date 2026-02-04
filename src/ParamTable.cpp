#include "ParamTable.h"

// Helper: all parameter names in pt that are NOT in premap_trend_params
std::unordered_set<std::string>
make_non_premap_param_set(const ParamTable& pt,
                          const std::unordered_set<std::string>& premap_trend_params)
{
  std::unordered_set<std::string> out;
  const int p = pt.base_names.size();
  out.reserve(p);

  for (int j = 0; j < p; ++j) {
    std::string nm = Rcpp::as<std::string>(pt.base_names[j]);
    if (premap_trend_params.find(nm) == premap_trend_params.end()) {
      out.insert(nm);
    }
  }
  return out;
}

