#include <unordered_map>
#include <memory>
#include <array>
#include <vector>     //
#include <Rcpp.h>    //
// #include "EMC2/userfun.hpp"

// View
struct KernelParsView {
  int n_rows;
  std::vector<const double*> cols;  // cols[k][row] = value for param k at trial row
};

// ---- Types ----

enum class KernelType {
  SimpleDelta,
  Delta2Kernel,
  Delta2LR,
  LinIncr,
  LinDecr,
  ExpIncr,
  ExpDecr,
  PowIncr,
  PowDecr,
  Poly2,
  Poly3,
  Poly4
  // Custom
};

// Some meta-data for kernels -- mostly for the future
struct KernelMeta {
  int  input_arity;          // how many *inputs* the kernel expects at once
  bool supports_grouping;    // whether a vector of names should be expanded into separate kernels
};

inline KernelMeta kernel_meta(KernelType kt) {
  switch (kt) {
  case KernelType::SimpleDelta:
  case KernelType::Delta2Kernel:
  case KernelType::Delta2LR:
  case KernelType::LinIncr:
  case KernelType::LinDecr:
  case KernelType::ExpIncr:
  case KernelType::ExpDecr:
  case KernelType::PowIncr:
  case KernelType::PowDecr:
  case KernelType::Poly2:
  case KernelType::Poly3:
  case KernelType::Poly4:
    return {1, true};   // all current kernels: 1D input, grouping allowed
  }

  // default future behaviour: 1D, but no grouping
  return {1, false};
}

// ---- Base + hierarchy ----
struct BaseKernel {
protected:
  std::vector<double> out_;
  bool has_run_ = false;

public:
  virtual ~BaseKernel() {}

  virtual void run(const KernelParsView& kernel_pars,
                   const Rcpp::NumericVector& covariate,
                   const std::vector<int>& comp_idx) = 0;

  virtual void reset() {
    out_.clear();
    has_run_ = false;
  }

  bool has_run() const { return has_run_; }

  const std::vector<double>& get_output() const { return out_; }

  // Expand compressed out_ (length n_comp) into full length using expand_idx
  void do_expand(const std::vector<int>& expand_idx) {
    const int n_full = static_cast<int>(expand_idx.size());

    out_.resize(n_full);  // keeps compressed data in [0..n_comp-1]

    for (int i = n_full - 1; i >= 0; --i) {
      int k = expand_idx[i] - 1;  // 1-based -> 0-based
      out_[i] = out_[k];
    }
  }

protected:
  void mark_run_complete() { has_run_ = true; }
};

// struct CustomKernel {
// protected:
//   std::vector<double> out_;
//   bool has_run_ = false;
//
// private:
//   Rcpp::XPtr<userfun_t> fun_;
//
// public:
//   CustomKernel(SEXP funptrSEXP)
//     : fun_(funptrSEXP) {
//     if (fun_.get() == nullptr) {
//       Rcpp::stop("CustomKernel: null function pointer.");
//     }
//     if (!(*fun_)) {
//       Rcpp::stop("CustomKernel: invalid function pointer.");
//     }
//   }
//
//   void run(const KernelParsView& kernel_pars,
//            const Rcpp::NumericVector& covariate,
//            const KernelParsView& par_input) {
//
//              // Call user function
//              userfun_t f = *fun_;
//              out_ = f(kernel_pars, covariate);
//
//              mark_run_complete();
//            }
// };

// For sequential kernels: currently same as BaseKernel -- just included to allow for other types (e.g. Bayesian ideal observer, autoregressive) in the future
struct SequentialKernel : BaseKernel {
  virtual ~SequentialKernel() {}
};

// All 1D delta kernels have scalar q and 1D pes_
struct DeltaKernel : SequentialKernel {
protected:
  double q_ = NA_REAL;             // latest value
  double pe_ = NA_REAL;            // latest PE
  std::vector<double> pes_;        // PE per trial

public:
  virtual ~DeltaKernel() {}

  const std::vector<double>& get_pes() const {
    return pes_;
  }
};


// ---- Individual kernels ----
// ---- Non-sequential kernels ----

struct LinIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);    // compressed output

             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 out_[j] = x;  // compressed index
               }
             }

             mark_run_complete();
           }
};


struct LinDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 out_[j] = -x;
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct ExpDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("ExpDecrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* lambda_col = kernel_pars.cols[0];
             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double lambda = lambda_col[r];
                 out_[j] = std::exp(-lambda * x);
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct ExpIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("ExpIncrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* lambda_col = kernel_pars.cols[0];
             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double lambda = lambda_col[r];
                 out_[j] = 1.0 - std::exp(-lambda * x);
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct PowDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("PowDecrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* alpha_col = kernel_pars.cols[0];
             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double alpha = alpha_col[r];
                 out_[j] = std::pow(1.0 + x, -alpha);
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct PowIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("PowIncrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* alpha_col = kernel_pars.cols[0];
             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double alpha = alpha_col[r];
                 out_[j] = 1.0 - std::pow(1.0 + x, -alpha);
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct Poly2Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("Poly2Kernel expects 2 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];

             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double a1 = a1_col[r];
                 double a2 = a2_col[r];
                 double x2 = x * x;
                 out_[j] = a1 * x + a2 * x2;
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct Poly3Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("Poly3Kernel expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];
             const double* a3_col = kernel_pars.cols[2];

             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double a1 = a1_col[r];
                 double a2 = a2_col[r];
                 double a3 = a3_col[r];
                 double x2 = x * x;
                 double x3 = x2 * x;
                 out_[j] = a1 * x + a2 * x2 + a3 * x3;
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct Poly4Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 4) {
               Rcpp::stop("Poly4Kernel expects 4 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];
             const double* a3_col = kernel_pars.cols[2];
             const double* a4_col = kernel_pars.cols[3];

             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double a1 = a1_col[r];
                 double a2 = a2_col[r];
                 double a3 = a3_col[r];
                 double a4 = a4_col[r];

                 double x2 = x * x;
                 double x3 = x2 * x;
                 double x4 = x2 * x2;
                 out_[j] = a1 * x + a2 * x2 + a3 * x3 + a4 * x4;
               }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};


// Sequential kernels
struct SimpleDelta : DeltaKernel {
  SimpleDelta() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("SimpleDelta expects 2 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             using namespace Rcpp;

             const int n_comp = static_cast<int>(comp_idx.size());
             if (n_comp <= 0) {
               out_.clear();
               pes_.clear();
               mark_run_complete();
               return;
             }

             out_.assign(n_comp, NA_REAL);
             pes_.assign(n_comp, NA_REAL);

             const double* q0_col    = kernel_pars.cols[0];
             const double* alpha_col = kernel_pars.cols[1];

             // Initial compressed element
             int row0 = comp_idx[0];             // full-data row index
             if (row0 < 0 || row0 >= covariate.size()) {
               stop("SimpleDelta::run: comp_idx[0] = %d out of range [0,%d)",
                    row0, covariate.size());
             }

             q_      = q0_col[row0];
             out_[0] = q_;

             double pe = NA_REAL;

             // j runs over COMPRESSED indices: 0..n_comp-2
             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];            // full-data row index
               if (r < 0 || r >= covariate.size()) {
                 stop("SimpleDelta::run: comp_idx[%d] = %d out of range [0,%d)",
                      j, r, covariate.size());
               }

               double x = covariate[r];
               if (!ISNAN(x)) {
                 double alpha = alpha_col[r];
                 pe = x - q_;
                 q_ += alpha * pe;
               } else {
                 pe = NA_REAL;
               }
               pes_[j] = pe;                   // compressed index

               int next_j = j + 1;
               out_[next_j] = q_;
             }

             mark_run_complete();
           }
};

struct Delta2LR : DeltaKernel {
  Delta2LR() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("Delta2LR expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, NA_REAL);
             pes_.assign(n_comp, NA_REAL);

             const double* q0_col       = kernel_pars.cols[0];
             const double* alphaPos_col = kernel_pars.cols[1];
             const double* alphaNeg_col = kernel_pars.cols[2];

             int row0 = comp_idx[0];
             out_[0] = q_ = q0_col[row0];

             double pe = NA_REAL;

             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               if (!ISNAN(x)) {
                 double alphaPos = alphaPos_col[r];
                 double alphaNeg = alphaNeg_col[r];
                 pe = x - q_;
                 double alpha = (pe > 0.0) ? alphaPos : alphaNeg;
                 q_ += alpha * pe;
               } else {
                 pe = NA_REAL;
               }

               pes_[j] = pe;           // compressed index
               out_[j + 1] = q_;
             }

             mark_run_complete();
           }
};

// 2D PE kernel: separate from DeltaKernel
struct Delta2Kernel : SequentialKernel {
  double qFast_ = NA_REAL;
  double qSlow_ = NA_REAL;
  double q_     = NA_REAL;

  // [compressed trial][0 = fast PE, 1 = slow PE]
  std::vector<std::array<double, 2>> pes2_;

  Delta2Kernel() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 4) {
               Rcpp::stop("Delta2Kernel expects 4 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, NA_REAL);
             pes2_.assign(n_comp, std::array<double,2>{NA_REAL, NA_REAL});

             const double* q0_col        = kernel_pars.cols[0];
             const double* alphaFast_col = kernel_pars.cols[1];
             const double* propSlow_col  = kernel_pars.cols[2];
             const double* dSwitch_col   = kernel_pars.cols[3];

             int row0 = comp_idx[0];
             out_[0] = qFast_ = qSlow_ = q_ = q0_col[row0];

             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];
               double x = covariate[r];
               double peFast = NA_REAL;
               double peSlow = NA_REAL;

               if (!ISNAN(x)) {
                 double alphaFast = alphaFast_col[r];
                 double propSlow  = propSlow_col[r];
                 double dSwitch   = dSwitch_col[r];
                 double alphaSlow = propSlow * alphaFast;

                 peFast = x - qFast_;
                 peSlow = x - qSlow_;

                 qFast_ += alphaFast * peFast;
                 qSlow_ += alphaSlow * peSlow;

                 double diff = std::abs(qFast_ - qSlow_);
                 q_ = (diff > dSwitch) ? qFast_ : qSlow_;
               }

               pes2_[j][0] = peFast;  // compressed index
               pes2_[j][1] = peSlow;
               out_[j + 1] = q_;
             }

             mark_run_complete();
           }
};

// ---- Type mapping + factory ----

KernelType to_kernel_type(const Rcpp::String& k);

std::unique_ptr<BaseKernel> make_kernel(KernelType kt);



