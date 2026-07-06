#include "kernels.h"

KernelType to_kernel_type(const Rcpp::String& k) {
  if (k == "delta")        return KernelType::SimpleDelta;
  if (k == "delta2lr")     return KernelType::Delta2LR;
  if (k == "delta_decoupled")     return KernelType::DeltaDecoupled;
  if (k == "delta2kernel") return KernelType::Delta2Kernel;
  if (k == "lin_incr")     return KernelType::LinIncr;
  if (k == "lin_decr")     return KernelType::LinDecr;
  if (k == "exp_incr")     return KernelType::ExpIncr;
  if (k == "exp_decr")     return KernelType::ExpDecr;
  if (k == "pow_incr")     return KernelType::PowIncr;
  if (k == "pow_decr")     return KernelType::PowDecr;
  if (k == "poly2")        return KernelType::Poly2;
  if (k == "poly3")        return KernelType::Poly3;
  if (k == "poly4")        return KernelType::Poly4;
  if (k == "custom")       return KernelType::Custom;
  if (k == "rescorlawagner")      return KernelType::RescorlaWagner;
  if (k == "beta_binomial")       return KernelType::BetaBinomial;
  if (k == "beta_binomial_decay") return KernelType::BetaBinomialDecay;
  if (k == "beta_binomial_window")return KernelType::BetaBinomialWindow;
  if (k == "dbm")                 return KernelType::DBM;
  if (k == "tpm")                 return KernelType::TPM;

  Rcpp::stop("Unknown kernel type");
}


std::unique_ptr<BaseKernel> make_kernel(KernelType kt, SEXP custom_fun) {
  switch (kt) {
  case KernelType::SimpleDelta:    return std::make_unique<SimpleDelta>();
  case KernelType::Delta2Kernel:   return std::make_unique<Delta2Kernel>();
  case KernelType::Delta2LR:       return std::make_unique<Delta2LR>();
  case KernelType::DeltaDecoupled: return std::make_unique<DeltaDecoupled>();

  case KernelType::LinIncr:     return std::make_unique<LinIncrKernel>();
  case KernelType::LinDecr:     return std::make_unique<LinDecrKernel>();
  case KernelType::ExpIncr:     return std::make_unique<ExpIncrKernel>();
  case KernelType::ExpDecr:     return std::make_unique<ExpDecrKernel>();
  case KernelType::PowIncr:     return std::make_unique<PowIncrKernel>();
  case KernelType::PowDecr:     return std::make_unique<PowDecrKernel>();
  case KernelType::Poly2:       return std::make_unique<Poly2Kernel>();
  case KernelType::Poly3:       return std::make_unique<Poly3Kernel>();
  case KernelType::Poly4:       return std::make_unique<Poly4Kernel>();
  case KernelType::RescorlaWagner:     return std::make_unique<RescorlaWagnerKernel>();
  case KernelType::BetaBinomial:       return std::make_unique<BetaBinomialKernel>();
  case KernelType::BetaBinomialDecay:  return std::make_unique<BetaBinomialDecayKernel>();
  case KernelType::BetaBinomialWindow: return std::make_unique<BetaBinomialWindowKernel>();
  case KernelType::DBM:                return std::make_unique<DBMKernel>();
  case KernelType::TPM:                return std::make_unique<TPMKernel>();

  case KernelType::Custom:
    if (custom_fun == R_NilValue) {
      Rcpp::stop("make_kernel: Custom kernel requested but custom_fun is NULL");
    }
    return std::unique_ptr<BaseKernel>(new CustomKernel(custom_fun));
  }

  Rcpp::stop("Unknown kernel type");
}
