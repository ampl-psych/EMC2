#ifndef EMC2_LIKELIHOOD_CONTEXT_H
#define EMC2_LIKELIHOOD_CONTEXT_H

#include <Rcpp.h>
#include <memory>

namespace emc2 {

// Copies share fixed inputs and own their scratch. Construct and destroy on
// the R thread; evaluation of thread-safe contexts uses only native memory.
class LikelihoodContext {
public:
  LikelihoodContext(Rcpp::CharacterVector parameter_names,
                    Rcpp::DataFrame data, Rcpp::NumericVector constants,
                    Rcpp::List designs, Rcpp::String type, Rcpp::List bounds,
                    Rcpp::List transforms, Rcpp::List pretransforms,
                    double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                    int n_threads = 1);
  ~LikelihoodContext();
  LikelihoodContext(LikelihoodContext&&) noexcept;
  LikelihoodContext& operator=(LikelihoodContext&&) noexcept;
  LikelihoodContext(const LikelihoodContext&);
  LikelihoodContext& operator=(const LikelihoodContext&) = delete;

  // Column-major coefficients. Output has one row, or trial_count() rows.
  // Parameter order is fixed by the constructor; n_particles may vary.
  void evaluate(const double* particles, int n_particles, double* result,
                bool trialwise = false);
  int trial_count() const;
  bool thread_safe() const;
  std::string unsafe_reason() const;   // why not, or ""

private:
  struct Impl;
  struct Workspace;
  std::shared_ptr<const Impl> impl_;
  std::unique_ptr<Workspace> work_;
};

} // namespace emc2
#endif
