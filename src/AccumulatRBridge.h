#ifndef EMC2_ACCUMULATR_BRIDGE_H
#define EMC2_ACCUMULATR_BRIDGE_H

#include <AccumulatR/api.hpp>
#include <algorithm>
#include <memory>
#include "ParamTable.h"
#include "transform_utils.h"

namespace emc2 {

// Fixed bindings and observation plans are shared; each chain owns its workspace.
struct AccumulatRBridge {
  struct Binding { int destination, source, length; };
  const accumulatr::LikelihoodAPI& api;
  std::unique_ptr<accumulatr::PreparedLikelihood,
                  void (*)(accumulatr::PreparedLikelihood*)> prepared;
  std::vector<Binding> bindings;
  std::vector<double> defaults;
  std::vector<int> counts;
  int rows_per_trial;

  struct Workspace {
    std::unique_ptr<accumulatr::LikelihoodWorkspace,
                    void (*)(accumulatr::LikelihoodWorkspace*)> native;
    std::vector<double> parameters;
    std::vector<int> valid;
    explicit Workspace(const AccumulatRBridge& bridge)
      : native(bridge.api.workspace(bridge.prepared.get()), bridge.api.destroy_workspace),
        parameters(bridge.defaults), valid(bridge.counts.size()) {}
  };

  AccumulatRBridge(Rcpp::DataFrame data, const ParamTable& pt,
                   std::vector<BoundSpec>& bounds, double min_ll)
    : api(accumulatr::likelihood_api()), prepared(nullptr, api.destroy) {
    Rcpp::List context = data.attr("AccumulatR_bridge");
    SEXP native = context["native"];
    if (TYPEOF(native) != EXTPTRSXP || !R_ExternalPtrAddr(native))
      Rcpp::stop("AccumulatR likelihood context is missing");
    Rcpp::List recipe = context["bridge"];
    Rcpp::NumericMatrix values = recipe["defaults"];
    Rcpp::CharacterMatrix sources = recipe["source_names"];
    defaults.assign(values.begin(), values.end());
    counts = Rcpp::as<std::vector<int>>(context["trial_counts"]);
    rows_per_trial = values.nrow() / counts.size();
    std::vector<unsigned char> used(pt.base.data.size(), 0);
    for (int col = 0; col < values.ncol(); ++col) {
      for (int row = 0; row < values.nrow(); ++row) {
        if (sources(row, col) == NA_STRING) continue;
        const int destination = row + values.nrow() * col;
        const int source = row + values.nrow() *
          pt.base_index_for(Rcpp::as<std::string>(sources(row, col)));
        used[source] = 1;
        if (!bindings.empty() && bindings.back().destination + bindings.back().length == destination &&
            bindings.back().source + bindings.back().length == source) ++bindings.back().length;
        else bindings.push_back({destination, source, 1});
      }
    }
    // Bounds only apply where a parameter is consumed by an active accumulator.
    std::vector<BoundSpec> active;
    for (const auto& bound : bounds) {
      const int offset = values.nrow() * bound.col_idx;
      int row = 0;
      while (row < values.nrow()) {
        if (!used[offset + row]) { ++row; continue; }
        BoundSpec block = bound;
        block.begin = row;
        while (row < values.nrow() && used[offset + row]) ++row;
        block.end = row;
        active.push_back(block);
      }
    }
    bounds = std::move(active);
    prepared.reset(api.prepare(native, data, values.nrow(), values.ncol(), min_ll));
  }

  void evaluate(const ParamTable& pt, const std::vector<int>& row_ok,
                 Workspace& work, double* out) const {
    for (const auto& binding : bindings)
      std::copy_n(pt.base.data.data() + binding.source, binding.length,
                  work.parameters.data() + binding.destination);
    bool all_valid = true;
    for (size_t trial = 0; trial < counts.size(); ++trial) {
      int valid = 1;
      for (int row = 0; row < rows_per_trial; ++row)
        valid &= row_ok[trial * rows_per_trial + row];
      work.valid[trial] = valid;
      all_valid = all_valid && valid;
    }
    api.evaluate(prepared.get(), work.native.get(), work.parameters.data(),
                  all_valid ? nullptr : work.valid.data(), out);
  }
};

} // namespace emc2
#endif
