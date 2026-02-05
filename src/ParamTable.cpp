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

// Helper: all parameter names in pt that are NOT in premap_trend_params OR in pretransform_trend_params
std::unordered_set<std::string>
  make_non_premap_pretransform_param_set(const ParamTable& pt,
                                         const std::unordered_set<std::string>& premap_trend_params,
                                         const std::unordered_set<std::string>& pretransform_trend_params)
  {
    std::unordered_set<std::string> out;
    const int p = pt.base_names.size();
    out.reserve(p);

    for (int j = 0; j < p; ++j) {
      std::string nm = Rcpp::as<std::string>(pt.base_names[j]);
      // skip if in *either* premap or pretransform set
      if (premap_trend_params.find(nm) != premap_trend_params.end()) continue;
      if (pretransform_trend_params.find(nm) != pretransform_trend_params.end()) continue;
      out.insert(std::move(nm));
    }
    return out;
  }

