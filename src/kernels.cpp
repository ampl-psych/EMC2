#include "kernels.h"

KernelType to_kernel_type(const Rcpp::String& k) {
  if (k == "delta2kernel") return KernelType::Delta2Kernel;
  if (k == "delta")        return KernelType::SimpleDelta;
  if (k == "delta2lr")     return KernelType::Delta2LR;
  if (k == "lin_incr")     return KernelType::LinIncr;
  if (k == "lin_decr")     return KernelType::LinDecr;
  if (k == "exp_incr")     return KernelType::ExpIncr;
  if (k == "exp_decr")     return KernelType::ExpDecr;
  if (k == "pow_incr")     return KernelType::PowIncr;
  if (k == "pow_decr")     return KernelType::PowDecr;
  if (k == "poly2")        return KernelType::Poly2;
  if (k == "poly3")        return KernelType::Poly3;
  if (k == "poly4")        return KernelType::Poly4;

  Rcpp::stop("Unknown kernel type");
}


std::unique_ptr<BaseKernel> make_kernel(KernelType kt) {
  switch (kt) {
  case KernelType::SimpleDelta: return std::unique_ptr<BaseKernel>(new SimpleDelta());
  case KernelType::Delta2Kernel:return std::unique_ptr<BaseKernel>(new Delta2Kernel());
  case KernelType::Delta2LR:    return std::unique_ptr<BaseKernel>(new Delta2LR());

  case KernelType::LinIncr:     return std::unique_ptr<BaseKernel>(new LinIncrKernel());
  case KernelType::LinDecr:     return std::unique_ptr<BaseKernel>(new LinDecrKernel());
  case KernelType::ExpIncr:     return std::unique_ptr<BaseKernel>(new ExpIncrKernel());
  case KernelType::ExpDecr:     return std::unique_ptr<BaseKernel>(new ExpDecrKernel());
  case KernelType::PowIncr:     return std::unique_ptr<BaseKernel>(new PowIncrKernel());
  case KernelType::PowDecr:     return std::unique_ptr<BaseKernel>(new PowDecrKernel());
  case KernelType::Poly2:       return std::unique_ptr<BaseKernel>(new Poly2Kernel());
  case KernelType::Poly3:       return std::unique_ptr<BaseKernel>(new Poly3Kernel());
  case KernelType::Poly4:       return std::unique_ptr<BaseKernel>(new Poly4Kernel());
  }
  Rcpp::stop("Unknown kernel type");
}
