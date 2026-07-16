#ifndef KERNELS_H
#define KERNELS_H

#include <unordered_map>
#include <memory>
#include <array>
#include <vector>     //
#include <Rcpp.h>    //
#include "nan_check.h"
#include "EMC2/userfun.hpp"
#include "Mat.h"
#include "kernels_math.h"

// View
struct KernelParsView {
  int n_rows;
  std::vector<const double*> cols;  // cols[k][row] = value for param k at trial row
};

// Struct for optional kernel arguments.
struct KernelArgs {
  const int* q_reset = nullptr;  // raw pointer into an IntegerVector; null = no reset
  int grid_res = 100; // resolution of discretised probability dists for DBM-like kernels
  const int* belief_reset = nullptr;
  // Future extensible fields go here, e.g.:
  // const double* some_other_col = nullptr;
};

// Struct for outputs. Outputs sometimes have 1 column, sometimes multiple. Pure C++ (for future threadsafe ops)
struct KernelOutput {
  const double* data = nullptr;  // raw pointer into kernel-owned storage
  int n_rows = 0;
  int n_cols = 1;                // 1 for all current kernels, N for RescorlaWagner

  // Convenience: element access (column-major, matches R matrix layout)
  double operator()(int r, int col) const {
    return data[col * n_rows + r];
  }
};

// ---- Types ----

enum class KernelType {
  SimpleDelta,
  Delta2Kernel,
  DeltaDecoupled,
  // Delta2Kernel2,
  Delta2LR,
  LinIncr,
  LinDecr,
  ExpIncr,
  ExpDecr,
  PowIncr,
  PowDecr,
  Poly2,
  Poly3,
  Poly4,
  Custom,
  RescorlaWagner,
  BetaBinomial,
  BetaBinomialDecay,
  BetaBinomialWindow,
  DBM,
  TPM
};

// Some meta-data for kernels -- mostly for the future
struct KernelMeta {
  int  input_arity;          // how many *inputs* the kernel expects at once
  bool supports_grouping;    // whether a vector of names should be expanded into separate kernels
};

inline KernelMeta kernel_meta(KernelType kt) {
  switch (kt) {
  case KernelType::SimpleDelta:
  case KernelType::DeltaDecoupled:
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
    return {1, true};   // all above kernels: 1D input, grouping allowed
  case KernelType::Custom: return{1, false};
  case KernelType::RescorlaWagner: return{-1, false};  // N columns allowed
  case KernelType::BetaBinomial:
  case KernelType::BetaBinomialDecay:
  case KernelType::BetaBinomialWindow:
  case KernelType::DBM:
  case KernelType::TPM:
    return {1, false};
  }

  // default future behaviour: 1D, but no grouping
  return {1, false};
}

// ---- Base + hierarchy ----
struct BaseKernel {
protected:
  std::vector<double> out_;
  bool has_run_ = false;

  // Remember expansion mapping (for 'at')
  std::vector<int> expand_idx_;   // 1-based indices
  bool has_expand_idx_ = false;

  // Two separate transpose buffers, one per stream family.
  // stream_buf_[0] is used for stream code 1 (Q / primary output)
  // stream_buf_[1] is used for stream code 2 (PE / secondary output)
  // Subclasses that need more streams can add further buffers explicitly.
  mutable std::vector<double> stream_buf_[2];

public:
  virtual ~BaseKernel() {}

  virtual void set_kernel_args(const KernelArgs& /*args*/) {}

  virtual void run(const KernelParsView& kernel_pars,
                   const Mat& covariate,
                   const std::vector<int>& comp_idx) = 0;

  virtual void reset() {
    out_.clear();
    stream_buf_[0].clear();
    stream_buf_[1].clear();
    has_run_ = false;
    expand_idx_.clear();
    has_expand_idx_ = false;
  }


  bool has_run() const { return has_run_; }

  // const std::vector<double>& get_output() const { return out_; }

  void set_expand_idx(const std::vector<int>& idx) {
    expand_idx_ = idx;
    has_expand_idx_ = !expand_idx_.empty();
  }
  const std::vector<int>& expand_idx() const { return expand_idx_; }
  bool has_expand_idx() const { return has_expand_idx_; }

  // Expand compressed out_ (length n_comp) into full length using expand_idx
  void do_expand(const std::vector<int>& expand_idx) {
    const int n_full = static_cast<int>(expand_idx.size());

    out_.resize(n_full);  // keeps compressed data in [0..n_comp-1]

    for (int i = n_full - 1; i >= 0; --i) {
      int k = expand_idx[i] - 1;  // 1-based -> 0-based
      out_[i] = out_[k];
    }
  }

  // does this kernel have a stream for this code?
  virtual bool has_output_stream(int code) const {
    return (code == 1);  // default: only main trajectory
  }

  // Returns a KernelOutput view into kernel-owned storage.
  // For single-column kernels: points directly into out_ (zero-copy).
  // For multi-column kernels (RescorlaWagner): points into col_major_buf_
  // which is populated lazily on first call.
  virtual KernelOutput get_output_stream(int code) const {
    if (code != 1) {
      Rcpp::stop("BaseKernel::get_output_stream: unsupported code %d (only 1)", code);
    }
    KernelOutput ko;
    ko.data   = out_.data();
    ko.n_rows = static_cast<int>(out_.size());
    ko.n_cols = 1;
    return ko;
  }

  // // single-stream getter, code=1 for main trajectory by default. code=2 for pes in delta, code=3 for xx in new kernels
  // virtual Rcpp::NumericVector get_output_stream(int code) const {
  //   using namespace Rcpp;
  //   if (code != 1) {
  //     stop("BaseKernel::get_output_stream: unsupported code %d (only 1)", code);
  //   }
  //   // out_ is already full-length at this point
  //   return wrap(out_);  // copies to NumericVector
  // }

  // Optional: name for each stream
  virtual std::string output_stream_name(int code) const {
    if (code == 1) return "covariate";
    throw std::runtime_error("BaseKernel::output_stream_name: unsupported code");
  }


protected:
  void mark_run_complete() { has_run_ = true; }
};

struct CustomKernel : BaseKernel {
private:
  Rcpp::XPtr<userfun_t> fun_;

public:
  // funptrSEXP is the external pointer stored in trend$custom_ptr
  CustomKernel(SEXP funptrSEXP) : fun_(funptrSEXP) {
    if (fun_.get() == nullptr) {
      Rcpp::stop("CustomKernel: null function pointer.");
    }
    if (!(*fun_)) {
      Rcpp::stop("CustomKernel: invalid function pointer.");
    }
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& input,
           const std::vector<int>& comp_idx) override {

             const int n_comp   = static_cast<int>(comp_idx.size());
             const int n_pars   = static_cast<int>(kernel_pars.cols.size());
             const int n_inputs = input.ncol;

             if (n_comp == 0) {
               out_.clear();
               mark_run_complete();
               return;
             }

             // 1) Build compressed parameter matrix: n_comp x n_pars
             Rcpp::NumericMatrix pars_comp(n_comp, n_pars);
             for (int p = 0; p < n_pars; ++p) {
               const double* col = kernel_pars.cols[p];
               for (int j = 0; j < n_comp; ++j) {
                 int r = comp_idx[j];        // full trial index
                 pars_comp(j, p) = col[r];   // compressed row j
               }
             }

             // 2) Build compressed input matrix: n_comp x n_inputs
             Rcpp::NumericMatrix input_comp(n_comp, n_inputs);
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               for (int c = 0; c < n_inputs; ++c) {
                 input_comp(j, c) = input(r, c);
               }
             }

             // 3) Call user function
             userfun_t f = *fun_;
             Rcpp::NumericVector res = f(pars_comp, input_comp);

             if (res.size() != n_comp) {
               Rcpp::stop("CustomKernel: user function returned length %d, expected %d (compressed trials)",
                          res.size(), n_comp);
             }

             // 4) Copy into out_ (compressed trajectory)
             out_.assign(n_comp, 0.0);
             for (int j = 0; j < n_comp; ++j) {
               out_[j] = res[j];
             }

             mark_run_complete();
           }
};

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
  const int* q_reset_ = nullptr;   // <-- ADD: null = no reset

public:
  virtual ~DeltaKernel() {}

  void set_kernel_args(const KernelArgs& args) override {
    q_reset_ = args.q_reset;
  }

  // const std::vector<double>& get_pes() const {
  //   return pes_;
  // }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 2);
  }

  KernelOutput get_output_stream(int code) const override {
    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      return KernelOutput{ out_.data(), n_full, 1 };
    }

    if (code == 2) {
      stream_buf_[1].resize(n_full);
      if (!has_expand_idx_) {
        if ((int)pes_.size() != n_full)
          Rcpp::stop("DeltaKernel: pes_ length mismatch");
        for (int i = 0; i < n_full; ++i) stream_buf_[1][i] = pes_[i];
      } else {
        const auto& idx = expand_idx_;
        for (int i = 0; i < n_full; ++i) stream_buf_[1][i] = pes_[idx[i] - 1];
      }
      return KernelOutput{ stream_buf_[1].data(), n_full, 1 };
    }

    Rcpp::stop("DeltaKernel::get_output_stream: unsupported code %d (1=Q,2=PE)", code);
  }


  // Rcpp::NumericVector get_output_stream(int code) const override {
  //   using namespace Rcpp;
  //
  //   const int n_full = static_cast<int>(out_.size());
  //
  //   if (code == 1) {
  //     // main trajectory: already full-length
  //     return wrap(out_);
  //   }
  //
  //   if (code == 2) {
  //     NumericVector res(n_full);
  //
  //     if (!has_expand_idx_) {
  //       // no 'at': one-to-one
  //       if ((int)pes_.size() != n_full)
  //         stop("DeltaKernel: pes_ length mismatch");
  //       for (int i = 0; i < n_full; ++i) res[i] = pes_[i];
  //     } else {
  //       // with 'at': expand from compressed index
  //       const auto& idx = expand_idx_;
  //       if ((int)idx.size() != n_full)
  //         stop("DeltaKernel: expand_idx length mismatch");
  //       for (int i = 0; i < n_full; ++i) {
  //         int k = idx[i] - 1;  // compressed index
  //         res[i] = pes_[k];
  //       }
  //     }
  //     return res;
  //   }
  //
  //   stop("DeltaKernel::get_output_stream: unsupported code %d (1=Q,2=PE)", code);
  // }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qvalue";
    if (code == 2) return "PE";
    throw std::runtime_error("DeltaKernel::output_stream_name: unsupported code");
  }
};


// ---- Individual kernels ----
// ---- Non-sequential kernels ----

struct LinIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);    // compressed output

             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate(r,0);
               out_[j] = x;
               // if (!is_nan(x)) {
               //   out_[j] = x;  // compressed index
               // }
             }

             mark_run_complete();
           }
};


struct LinDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate(r,0);
               out_[j] = -x;
               // if (!is_nan(x)) {
               //   out_[j] = -x;
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct ExpDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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
               double x = covariate(r,0);
               // if (!is_nan(x)) {
                 double lambda = lambda_col[r];
                 out_[j] = std::exp(-lambda * x);
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct ExpIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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
               double x = covariate(r,0);
               // if (!is_nan(x)) {
                 double lambda = lambda_col[r];
                 out_[j] = 1.0 - std::exp(-lambda * x);
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct PowDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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
               double x = covariate(r,0);
               // if (!is_nan(x)) {
                 double alpha = alpha_col[r];
                 out_[j] = std::pow(1.0 + x, -alpha);
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct PowIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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
               double x = covariate(r,0);
               // if (!is_nan(x)) {
                 double alpha = alpha_col[r];
                 out_[j] = 1.0 - std::pow(1.0 + x, -alpha);
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct Poly2Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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
               double x = covariate(r,0);
               // if (!is_nan(x)) {
                 double a1 = a1_col[r];
                 double a2 = a2_col[r];
                 double x2 = x * x;
                 out_[j] = a1 * x + a2 * x2;
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct Poly3Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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
               double x = covariate(r,0);
               // if (!is_nan(x)) {
                 double a1 = a1_col[r];
                 double a2 = a2_col[r];
                 double a3 = a3_col[r];
                 double x2 = x * x;
                 double x3 = x2 * x;
                 out_[j] = a1 * x + a2 * x2 + a3 * x3;
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};

struct Poly4Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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
               double x = covariate(r,0);
               // if (!is_nan(x)) {
                 double a1 = a1_col[r];
                 double a2 = a2_col[r];
                 double a3 = a3_col[r];
                 double a4 = a4_col[r];

                 double x2 = x * x;
                 double x3 = x2 * x;
                 double x4 = x2 * x2;
                 out_[j] = a1 * x + a2 * x2 + a3 * x3 + a4 * x4;
               // }
               // out_[j] = last;
             }

             mark_run_complete();
           }
};


// Sequential kernels
struct SimpleDelta : DeltaKernel {
  SimpleDelta() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("SimpleDelta expects 2 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n_comp = static_cast<int>(comp_idx.size());
             if (n_comp <= 0) {
               out_.clear();
               pes_.clear();
               mark_run_complete();
               return;
             }

             // slightly faster -- no-op size check, no need to write anything
             out_.resize(n_comp);
             pes_.resize(n_comp);
             pes_[n_comp - 1] = NA_REAL;

             const double* q0_col    = kernel_pars.cols[0];
             const double* alpha_col = kernel_pars.cols[1];
             const double* cov_ptr   = covariate.colptr(0);

             int row0 = comp_idx[0];
             q_       = q0_col[row0];
             out_[0]  = q_;

             for (int j = 0; j < n_comp - 1; ++j) {
               int r    = comp_idx[j];
               // --- RESET (before PE) ---
               if (q_reset_ && q_reset_[r]) {
                 q_ = q0_col[r];
                 out_[j] = q_;           // overwrite with reset value
               }

               double x = cov_ptr[r];

               if (!is_nan(x)) {
                 double alpha = alpha_col[r];
                 double pe    = x - q_;
                 pes_[j]      = pe;
                 q_          += alpha * pe;
               } else {
                 pes_[j] = NA_REAL;
               }
               out_[j + 1] = q_;
             }

             mark_run_complete();
           }
};

// Delta rule reparametrised to decouple the movement towards the outcome from the decay towards 0
struct DeltaDecoupled : DeltaKernel {
  DeltaDecoupled() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("DeltaDecoupled expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n_comp = static_cast<int>(comp_idx.size());
             if (n_comp <= 0) {
               out_.clear();
               pes_.clear();
               mark_run_complete();
               return;
             }

             // slightly faster -- no-op size check, no need to write anything
             out_.resize(n_comp);
             pes_.resize(n_comp);
             pes_[n_comp - 1] = NA_REAL;

             const double* q0_col    = kernel_pars.cols[0];
             const double* alpha_col = kernel_pars.cols[1];
             const double* lambda_col = kernel_pars.cols[2];
             const double* cov_ptr   = covariate.colptr(0);

             int row0 = comp_idx[0];
             q_       = q0_col[row0];
             out_[0]  = q_;

             for (int j = 0; j < n_comp - 1; ++j) {
               int r    = comp_idx[j];
               // --- RESET (before PE) ---
               if (q_reset_ && q_reset_[r]) {
                 q_ = q0_col[r];
                 out_[j] = q_;           // overwrite with reset value
               }

               double x = cov_ptr[r];

               if (!is_nan(x)) {
                 double alpha  = alpha_col[r];
                 double lambda = lambda_col[r];
                 double term1 = alpha*x;
                 double term2 = lambda*q_;
                 q_ += term1 - term2;
                 double pe    = x - q_;  // bit questionable what the RPE is in this case though
                 pes_[j]      = pe;
                 // q_          += alpha * pe;
               } else {
                 pes_[j] = NA_REAL;
               }
               out_[j + 1] = q_;
             }

             mark_run_complete();
           }
};

struct Delta2LR : DeltaKernel {
  Delta2LR() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
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

               // --- RESET (before PE) ---
               if (q_reset_ && q_reset_[r]) {
                 q_ = q0_col[r];
                 out_[j] = q_;           // overwrite with reset value
               }

               double x = covariate(r,0);
               if (!is_nan(x)) {
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
  const int* q_reset_ = nullptr;

  void set_kernel_args(const KernelArgs& args) override {
    q_reset_ = args.q_reset;
  }

  // [compressed trial][0 = fast PE, 1 = slow PE]
  std::vector<double> pes_fast_;
  std::vector<double> pes_slow_;
  std::vector<double> q_fast_;
  std::vector<double> q_slow_;

  // One dedicated transpose buffer per secondary stream code.
  // Indexed as: secondary_buf_[code - 2], i.e.:
  //   code 2 (Qfast)   -> secondary_buf_[0]
  //   code 3 (Qslow)   -> secondary_buf_[1]
  //   code 4 (PEfast)  -> secondary_buf_[2]
  //   code 5 (PEslow)  -> secondary_buf_[3]
  mutable std::vector<double> secondary_buf_[4];

  Delta2Kernel() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 4) {
               Rcpp::stop("Delta2Kernel expects 4 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, NA_REAL);
             q_fast_.assign(n_comp, NA_REAL);
             q_slow_.assign(n_comp, NA_REAL);
             pes_fast_.assign(n_comp, NA_REAL);
             pes_slow_.assign(n_comp, NA_REAL);

             const double* q0_col        = kernel_pars.cols[0];
             const double* alphaFast_col = kernel_pars.cols[1];
             const double* propSlow_col  = kernel_pars.cols[2];
             const double* dSwitch_col   = kernel_pars.cols[3];

             int row0 = comp_idx[0];
             out_[0] = qFast_ = qSlow_ = q_ = q0_col[row0];

             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];

               // --- RESET (before PE): both trackers reset to q0 ---
               if (q_reset_ && q_reset_[r]) {
                 qFast_ = qSlow_ = q_ = q0_col[r];
                 out_[j] = q_;           // overwrite with reset value
               }

               double x = covariate(r,0);
               double peFast = NA_REAL;
               double peSlow = NA_REAL;

               if (!is_nan(x)) {
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
               q_fast_[j+1] = qFast_;  // compressed index
               q_slow_[j+1] = qSlow_;

               pes_fast_[j] = peFast;  // compressed index
               pes_slow_[j] = peSlow;
               out_[j + 1] = q_;
             }

             mark_run_complete();
           }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 5);
  }

  KernelOutput get_output_stream(int code) const override {
    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      return KernelOutput{ out_.data(), n_full, 1 };
    }

    const std::vector<double>* src = nullptr;
    if      (code == 2) src = &q_fast_;
    else if (code == 3) src = &q_slow_;
    else if (code == 4) src = &pes_fast_;
    else if (code == 5) src = &pes_slow_;
    else Rcpp::stop("Delta2Kernel::get_output_stream: unsupported code %d "
                      "(1=Q, 2=Qfast, 3=Qslow, 4=PEfast, 5=PEslow)", code);

    std::vector<double>& buf = secondary_buf_[code - 2];
    buf.resize(n_full);

    if (!has_expand_idx_) {
      if ((int)src->size() != n_full)
        Rcpp::stop("Delta2Kernel::get_output_stream: source length mismatch");
      for (int i = 0; i < n_full; ++i) buf[i] = (*src)[i];
    } else {
      const auto& idx = expand_idx_;
      for (int i = 0; i < n_full; ++i) buf[i] = (*src)[idx[i] - 1];
    }

    return KernelOutput{ buf.data(), n_full, 1 };
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qvalue";
    if (code == 2) return "Qfast";
    if (code == 3) return "Qslow";
    if (code == 4) return "PEfast";
    if (code == 5) return "PEslow";
    throw std::runtime_error("Delta2Kernel::output_stream_name: unsupported code");
  }
};


struct RescorlaWagnerKernel : SequentialKernel {
private:
  // Row-major internal storage: index as [r * n_covs_ + col]
  int n_covs_ = 0;
  std::vector<double> q_mat_;   // [n_comp * n_covs_]: Q-value per trial per covariate
  std::vector<double> pe_mat_;  // [n_comp * n_covs_]: compound PE for active covariates, NA otherwise

  const int* q_reset_ = nullptr;

public:
  void set_kernel_args(const KernelArgs& args) override {
    q_reset_ = args.q_reset;
  }

  void reset() override {
    BaseKernel::reset();
    q_mat_.clear();
    pe_mat_.clear();
    n_covs_ = 0;
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("RescorlaWagnerKernel expects 2 parameter columns (q0, alpha), got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n_comp = static_cast<int>(comp_idx.size());
             n_covs_          = covariate.ncol;

             if (n_comp == 0 || n_covs_ == 0) {
               q_mat_.clear();
               pe_mat_.clear();
               mark_run_complete();
               return;
             }

             const double* q0_col    = kernel_pars.cols[0];
             const double* alpha_col = kernel_pars.cols[1];

             q_mat_.assign(n_comp * n_covs_, NA_REAL);
             pe_mat_.assign(n_comp * n_covs_, NA_REAL);

             // Initialise: all covariates start at q0 of first trial
             int row0       = comp_idx[0];
             double q0_init = q0_col[row0];

             std::vector<double> q_cur(n_covs_, q0_init);

             // Write Q-values entering trial 0
             for (int c = 0; c < n_covs_; ++c) {
               q_mat_[0 * n_covs_ + c] = q_cur[c];
             }

             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];

               // Reset before PE: overwrite q_cur and the already-written q_mat_[j]
               if (q_reset_ && q_reset_[r]) {
                 double q0_r = q0_col[r];
                 for (int c = 0; c < n_covs_; ++c) {
                   q_cur[c]                  = q0_r;
                   q_mat_[j * n_covs_ + c]  = q0_r;  // overwrite entering Q for this trial
                 }
               }

               // Identify active covariates and accumulate compound Q
               double reward    = NA_REAL;
               double q_active  = 0.0;
               bool   any_active = false;

               for (int c = 0; c < n_covs_; ++c) {
                 double x = covariate(r, c);
                 if (!is_nan(x)) {
                   reward     = x;   // reward is the same across all active columns
                   q_active  += q_cur[c];
                   any_active = true;
                 }
               }

               if (any_active && !is_nan(reward)) {
                 double alpha       = alpha_col[r];
                 double compound_pe = reward - q_active;

                 for (int c = 0; c < n_covs_; ++c) {
                   double x = covariate(r, c);
                   if (!is_nan(x)) {
                     pe_mat_[j * n_covs_ + c] = compound_pe;
                     q_cur[c] += alpha * compound_pe;
                   }
                   // inactive: pe_mat_ stays NA, q_cur[c] unchanged
                 }
               }

               // Write Q-values entering trial j+1
               for (int c = 0; c < n_covs_; ++c) {
                 q_mat_[(j + 1) * n_covs_ + c] = q_cur[c];
               }
               // Rprintf("n_covs_=%d n_comp=%d\n", n_covs_, n_comp);
             }

             // pe_mat_ for the last trial stays NA (no outcome consumed yet),
             // mirroring SimpleDelta's pes_[n_comp - 1] = NA_REAL

             mark_run_complete();
           }

  bool has_output_stream(int code) const override {
    return (code == 1 || code == 2);
  }

  // Stream 1: Q-matrix (n_rows x n_covs_), column-major
  // Stream 2: PE-matrix (n_rows x n_covs_), column-major
  KernelOutput get_output_stream(int code) const override {
    if (code != 1 && code != 2) {
      Rcpp::stop("RescorlaWagnerKernel::get_output_stream: unsupported code %d "
                   "(1=Qmatrix, 2=PEmatrix)", code);
    }

    const std::vector<double>& src = (code == 1) ? q_mat_ : pe_mat_;
    const int n_comp  = static_cast<int>(src.size()) / n_covs_;
    std::vector<double>& buf = stream_buf_[code - 1];

    if (has_expand_idx_) {
      // Expand compact [n_comp x n_covs_] to full [n_trials x n_covs_],
      // then transpose to column-major for R.
      const int n_full = static_cast<int>(expand_idx_.size());
      buf.resize(n_full * n_covs_);

      for (int c = 0; c < n_covs_; ++c) {
        for (int i = 0; i < n_full; ++i) {
          int comp_row = expand_idx_[i] - 1;  // 1-based -> 0-based
          buf[c * n_full + i] = src[comp_row * n_covs_ + c];
        }
      }

      return KernelOutput{ buf.data(), n_full, n_covs_ };

    } else {
      // No expand: just transpose row-major [n_comp x n_covs_] to column-major
      buf.resize(n_comp * n_covs_);

      for (int c = 0; c < n_covs_; ++c) {
        for (int r = 0; r < n_comp; ++r) {
          buf[c * n_comp + r] = src[r * n_covs_ + c];
        }
      }

      return KernelOutput{ buf.data(), n_comp, n_covs_ };
    }
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qmatrix";
    if (code == 2) return "PEmatrix";
    throw std::runtime_error("RescorlaWagnerKernel::output_stream_name: unsupported code");
  }
};


// // 2kernel adjusted
// struct Delta2Kernel2 : Delta2Kernel {
//   double qFast_ = NA_REAL;
//   double qSlow_ = NA_REAL;
//   double q_     = NA_REAL;
//
//   Delta2Kernel2() {}
//
//   void run(const KernelParsView& kernel_pars,
//            const Rcpp::NumericMatrix& covariate,
//            const std::vector<int>& comp_idx) override {
//              if (kernel_pars.cols.size() != 4) {
//                Rcpp::stop("Delta2Kernel expects 4 parameter columns, got %d",
//                           (int)kernel_pars.cols.size());
//              }
//
//              int n_comp = comp_idx.size();
//              out_.assign(n_comp, NA_REAL);
//              q_fast_.assign(n_comp, NA_REAL);
//              q_slow_.assign(n_comp, NA_REAL);
//              pes_fast_.assign(n_comp, NA_REAL);
//              pes_slow_.assign(n_comp, NA_REAL);
//
//              const double* q0_col        = kernel_pars.cols[0];
//              const double* alphaFast_col = kernel_pars.cols[1];
//              const double* propSlow_col  = kernel_pars.cols[2];
//              const double* dSwitch_col   = kernel_pars.cols[3];
//
//              int row0 = comp_idx[0];
//              out_[0] = qFast_ = qSlow_ = q_ = q0_col[row0];
//              int current_kernel = 0; // 0 = fast, 1 = slow
//
//              for (int j = 0; j < n_comp - 1; ++j) {
//                int r = comp_idx[j];
//                double x = covariate(r,0);
//                double peFast = NA_REAL;
//                double peSlow = NA_REAL;
//
//                if (!ISNAN(x)) {
//                  double alphaFast = alphaFast_col[r];
//                  double propSlow  = propSlow_col[r];
//                  double dSwitch   = dSwitch_col[r];
//                  double alphaSlow = propSlow * alphaFast;
//
//                  peFast = x - qFast_;
//                  peSlow = x - qSlow_;
//
//                  qFast_ += alphaFast * peFast;
//                  qSlow_ += alphaSlow * peSlow;
//
//                  double diff = std::abs(qFast_ - qSlow_);
//                  if(diff > dSwitch) {
//                    current_kernel = 0; // fast kernel
//                    q_ = qFast_;
//                  } else {
//                    if(current_kernel == 0) {
//                      // was in fast mode, now moving to slow mode. Override Q-value of slow
//                      qSlow_ = qFast_;
//                    }
//                    current_kernel = 1;
//                    q_ = qSlow_;
//                  }
//                  // q_ = (diff > dSwitch) ? qFast_ : qSlow_;
//                }
//
//                q_fast_[j+1] = qFast_;  // compressed index
//                q_slow_[j+1] = qSlow_;
//
//                pes_fast_[j] = peFast;  // compressed index
//                pes_slow_[j] = peSlow;
//                out_[j + 1] = q_;
//              }
//
//              mark_run_complete();
//            }
// };

// =============================================================================
// DBMBaseKernel
// Streams: 1 = prediction mean, 2 = prediction mode, 3 = surprise (bits),
//          4 = prediction log-precision
// =============================================================================

struct DBMBaseKernel : BaseKernel {
protected:
  std::vector<double> pred_mean_;
  std::vector<double> pred_mode_;
  mutable std::vector<double> surprise_;          // computed lazily
  std::vector<double> pred_logprecision_;
  std::vector<double> comp_obs_;                  // compressed observations, stored during run()
  mutable bool surprise_computed_ = false;
  const int* belief_reset_ = nullptr;   // <-- ADD: null = no reset

  void store_obs(const double* cov_ptr, const std::vector<int>& comp_idx) {
    const int n_comp = static_cast<int>(comp_idx.size());
    comp_obs_.resize(n_comp);
    for (int j = 0; j < n_comp; ++j)
      comp_obs_[j] = cov_ptr[comp_idx[j]];
  }

  void ensure_surprise() const {
    if (surprise_computed_) return;
    const int n_comp = static_cast<int>(pred_mean_.size());
    surprise_.resize(n_comp, std::numeric_limits<double>::quiet_NaN());
    for (int j = 0; j < n_comp; ++j) {
      if (!is_nan(comp_obs_[j])) {
        surprise_[j] = shannon_surprise(pred_mean_[j], comp_obs_[j]);
      }
      // if (is_nan(comp_obs_[j])) {
      //   surprise_[j] = shannon_entropy(pred_mean_[j]);
      // } else {
      //   surprise_[j] = shannon_surprise(pred_mean_[j], comp_obs_[j]);
      // }
    }
    surprise_computed_ = true;
  }

public:
  void reset() override {
    BaseKernel::reset();
    pred_mean_.clear();
    pred_mode_.clear();
    surprise_.clear();
    pred_logprecision_.clear();
    comp_obs_.clear();
    surprise_computed_ = false;
  }

  void set_kernel_args(const KernelArgs& args) override {
    belief_reset_ = args.belief_reset;
  }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 4);
  }

  KernelOutput get_output_stream(int code) const override {
    if (code == 3) ensure_surprise(); // compute surprise the moment it is requested, not before
    const std::vector<double>* src = nullptr;
    if      (code == 1) src = &pred_mean_;
    else if (code == 2) src = &pred_mode_;
    else if (code == 3) src = &surprise_;
    else if (code == 4) src = &pred_logprecision_;
    else Rcpp::stop("DBMBaseKernel::get_output_stream: unsupported code %d "
                      "(1=mean, 2=mode, 3=surprise, 4=log-precision)", code);

    if (has_expand_idx_) {
      const int n_full = static_cast<int>(expand_idx_.size());
      stream_buf_[0].resize(n_full);
      for (int i = 0; i < n_full; ++i)
        stream_buf_[0][i] = (*src)[expand_idx_[i] - 1];
      return KernelOutput{ stream_buf_[0].data(), n_full, 1 };
    }
    return KernelOutput{ src->data(), static_cast<int>(src->size()), 1 };
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "mean";
    if (code == 2) return "mode";
    if (code == 3) return "surprise";
    if (code == 4) return "log-precision";
    throw std::runtime_error("DBMBaseKernel::output_stream_name: unsupported code");
  }

protected:
  void compute_surprise(const double* cov_ptr,
                        const std::vector<int>& comp_idx) {
    const int n_comp = static_cast<int>(comp_idx.size());
    surprise_.resize(n_comp, std::numeric_limits<double>::quiet_NaN());
    for (int j = 0; j < n_comp; ++j) {
      const double obs = cov_ptr[comp_idx[j]];
      if (!is_nan(obs)) {
        surprise_[j] = shannon_surprise(pred_mean_[j], obs);
      }
      // if (is_nan(obs)) {
      //   surprise_[j] = shannon_entropy(pred_mean_[j]);
      // } else {
      //   surprise_[j] = shannon_surprise(pred_mean_[j], obs);
      // }
    }
  }

  // keep out_ in sync with pred_mean_ so BaseKernel::do_expand works if called
  void sync_out_to_mean() { out_ = pred_mean_; }
};

// =============================================================================
// BetaBinomialKernel  —  basic (no memory constraint)
// Parameters: a0, b0
// =============================================================================

struct BetaBinomialKernel : DBMBaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 2)
               Rcpp::stop("BetaBinomialKernel expects 2 parameter columns (a0, b0), got %d",
                          (int)kernel_pars.cols.size());

             const int     n_comp  = static_cast<int>(comp_idx.size());
             const double* a0_col  = kernel_pars.cols[0];
             const double* b0_col  = kernel_pars.cols[1];
             const double* cov_ptr = covariate.colptr(0);

             pred_mean_.resize(n_comp);
             pred_mode_.resize(n_comp);
             pred_logprecision_.resize(n_comp);

             double n_hit = 0.0, n_trial = 0.0;

             for (int j = 0; j < n_comp; ++j) {
               const int    r   = comp_idx[j];

               if (belief_reset_ && belief_reset_[r]) {
                 double n_hit = 0.0, n_trial = 0.0;
               }

               const double a_t = a0_col[r] + n_hit;
               const double b_t = b0_col[r] + (n_trial - n_hit);

               pred_mean_[j] = beta_mean(a_t, b_t);
               pred_mode_[j] = beta_mode(a_t, b_t);
               pred_logprecision_[j] = beta_log_precision(a_t, b_t);

               const double x = cov_ptr[r];
               if (!is_nan(x)) { n_hit += x; n_trial += 1.0; }
             }

             store_obs(cov_ptr, comp_idx);
             sync_out_to_mean();
             mark_run_complete();
           }
};

// =============================================================================
// BetaBinomialDecayKernel  —  exponential decay on accumulated counts
// Parameters: a0, b0, decay
// =============================================================================

struct BetaBinomialDecayKernel : DBMBaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("BetaBinomialDecayKernel expects 3 parameter columns "
                            "(a0, b0, decay), got %d",
                            (int)kernel_pars.cols.size());

             const int     n_comp    = static_cast<int>(comp_idx.size());
             const double* a0_col    = kernel_pars.cols[0];
             const double* b0_col    = kernel_pars.cols[1];
             const double* decay_col = kernel_pars.cols[2];
             const double* cov_ptr   = covariate.colptr(0);

             pred_mean_.resize(n_comp);
             pred_mode_.resize(n_comp);
             pred_logprecision_.resize(n_comp);

             double n_hit = 0.0, n_trial = 0.0;

             for (int j = 0; j < n_comp; ++j) {
               const int    r   = comp_idx[j];

               if (belief_reset_ && belief_reset_[r]) {
                 double n_hit = 0.0, n_trial = 0.0;
               }

               const double a_t = a0_col[r] + n_hit;
               const double b_t = b0_col[r] + (n_trial - n_hit);

               pred_mean_[j] = beta_mean(a_t, b_t);
               pred_mode_[j] = beta_mode(a_t, b_t);
               pred_logprecision_[j] = beta_log_precision(a_t, b_t);

               const double df = std::exp(-1.0 / decay_col[r]);
               const double x  = cov_ptr[r];
               if (!is_nan(x)) {
                 n_hit   = df * (n_hit + x);
                 n_trial = df * (n_trial + 1.0);
               } else {
                 // still decay without observation: passage of time erodes memory
                 n_hit   = df * n_hit;
                 n_trial = df * n_trial;
               }
             }

             store_obs(cov_ptr, comp_idx);
             sync_out_to_mean();
             mark_run_complete();
           }
};

// =============================================================================
// BetaBinomialWindowKernel  —  fixed sliding window
// Parameters: a0, b0, window
// =============================================================================

struct BetaBinomialWindowKernel : DBMBaseKernel {
private:
  struct Event { double obs; int idx; };

public:
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("BetaBinomialWindowKernel expects 3 parameter columns "
                            "(a0, b0, window), got %d",
                            (int)kernel_pars.cols.size());

             const int     n_comp     = static_cast<int>(comp_idx.size());
             const double* a0_col     = kernel_pars.cols[0];
             const double* b0_col     = kernel_pars.cols[1];
             const double* window_col = kernel_pars.cols[2];
             const double* cov_ptr    = covariate.colptr(0);

             pred_mean_.resize(n_comp);
             pred_mode_.resize(n_comp);
             pred_logprecision_.resize(n_comp);

             double n_hit = 0.0, n_trial = 0.0;
             std::deque<Event> buf;

             for (int j = 0; j < n_comp; ++j) {
               const int r = comp_idx[j];
               const int w = static_cast<int>(window_col[r]);

               if (belief_reset_ && belief_reset_[r]) {
                 double n_hit = 0.0, n_trial = 0.0;
                 buf.clear();
               }

               // prune observations outside the window
               while (!buf.empty() && (r - buf.front().idx) > w) {
                 n_hit   -= buf.front().obs;
                 n_trial -= 1.0;
                 buf.pop_front();
               }

               const double a_t = a0_col[r] + n_hit;
               const double b_t = b0_col[r] + (n_trial - n_hit);

               pred_mean_[j] = beta_mean(a_t, b_t);
               pred_mode_[j] = beta_mode(a_t, b_t);
               pred_logprecision_[j] = beta_log_precision(a_t, b_t);

               const double x = cov_ptr[r];
               if (!is_nan(x)) {
                 buf.push_back({x, r});
                 n_hit   += x;
                 n_trial += 1.0;
               }
             }

             store_obs(cov_ptr, comp_idx);
             sync_out_to_mean();
             mark_run_complete();
           }
};

// =============================================================================
// DBMKernel  —  Dynamic Belief Model
// Yu & Cohen (2008), Ide et al. (2013)
// Parameters: cp, mu0, s0
// kernel_args: grid_res (default 100)
// =============================================================================

struct DBMKernel : DBMBaseKernel {
private:
  int grid_res_ = 100;

public:
  void set_kernel_args(const KernelArgs& args) override {
    if (args.grid_res > 0) grid_res_ = args.grid_res;
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("DBMKernel expects 3 parameter columns (cp, mu0, s0), got %d",
                          (int)kernel_pars.cols.size());

             const int     n_comp  = static_cast<int>(comp_idx.size());
             const double* cp_col  = kernel_pars.cols[0];
             const double* mu0_col = kernel_pars.cols[1];
             const double* s0_col  = kernel_pars.cols[2];
             const double* cov_ptr = covariate.colptr(0);

             pred_mean_.resize(n_comp);
             pred_mode_.resize(n_comp);
             pred_logprecision_.resize(n_comp);

             const int    gs     = grid_res_ + 1;

             std::vector<double> prob_grid(gs), x_like(gs), y_like(gs);
             for (int i = 0; i < gs; ++i) {
               prob_grid[i] = static_cast<double>(i) / (gs - 1);
               x_like[i]   = prob_grid[i];
               y_like[i]   = 1.0 - prob_grid[i];
             }

             std::vector<double> DBM_prior(gs), DBM_pred(gs), DBM_post(gs);

             for (int j = 0; j < n_comp; ++j) {
               const int    r   = comp_idx[j];
               const double cp  = cp_col[r];
               const double mu0 = mu0_col[r];
               const double s0  = s0_col[r];
               const double a   = mu0 * s0;
               const double b   = (1.0 - mu0) * s0;
               const double x   = cov_ptr[r];

               // compute discretised Beta prior
               for (int i = 0; i < gs; ++i)
                 DBM_prior[i] = dbeta_val(prob_grid[i], a, b);
               normalise_inplace(DBM_prior);

               // predictive distribution
               if (j == 0 || (belief_reset_ && belief_reset_[r])) {
                 // first trial or belief reset: fixed prior is the predictive
                 DBM_pred = DBM_prior;
               } else {
                 // otherwise: mixture of fixed prior and most recent posterior
                 const double mix_old = 1.0 - cp;
                 const double mix_new = cp;
                 for (int i = 0; i < gs; ++i)
                   DBM_pred[i] = mix_old * DBM_post[i] + mix_new * DBM_prior[i];
                 normalise_inplace(DBM_pred);
               }

               pred_mean_[j] = mean_discrete(prob_grid, DBM_pred);
               pred_mode_[j] = mode_discrete(prob_grid, DBM_pred);
               pred_logprecision_[j] = log_precision_discrete(prob_grid, DBM_pred);

               // posterior update
               if (is_nan(x)) {
                 DBM_post = DBM_pred;   // no observation: push predictive forward
               } else {
                 const std::vector<double>& like = (x == 1.0) ? x_like : y_like;
                 for (int i = 0; i < gs; ++i)
                   DBM_post[i] = DBM_pred[i] * like[i];
                 normalise_inplace(DBM_post);
               }
             }

             store_obs(cov_ptr, comp_idx);
             sync_out_to_mean();
             mark_run_complete();
           }
};

// =============================================================================
// TPMKernel  —  Transition Probability Model
// Meyniel et al. (2016)
// Parameters: cp, a0, b0
// kernel_args: grid_res (default 100)
// =============================================================================

struct TPMKernel : DBMBaseKernel {
private:
  int grid_res_ = 100;

  struct TPMGrid {
    int resol = 0, n_combi = 0;
    std::vector<double> p_XX, p_XY;
    std::vector<double> like_XX, like_XY, like_YX, like_YY;
    std::vector<double> mean_p;
  };

  TPMGrid build_grid(int grid_res) const {
    const int resol   = grid_res + 1;
    const int n_combi = resol * resol;

    std::vector<double> grid(resol);
    for (int i = 0; i < resol; ++i)
      grid[i] = static_cast<double>(i) / (resol - 1);

    TPMGrid g;
    g.resol = resol; g.n_combi = n_combi;
    g.p_XX.resize(n_combi);    g.p_XY.resize(n_combi);
    g.like_XX.resize(n_combi); g.like_XY.resize(n_combi);
    g.like_YX.resize(n_combi); g.like_YY.resize(n_combi);
    g.mean_p.resize(n_combi);

    int idx = 0;
    for (int i0 = 0; i0 < resol; ++i0) {
      const double pXY = grid[i0];
      for (int i1 = 0; i1 < resol; ++i1) {
        const double pXX   = grid[i1];
        g.p_XX[idx]    = pXX;
        g.p_XY[idx]    = pXY;
        g.like_XX[idx] = pXX;
        g.like_XY[idx] = pXY;
        g.like_YX[idx] = 1.0 - pXX;
        g.like_YY[idx] = 1.0 - pXY;
        g.mean_p[idx]  = 0.5 * (pXX + pXY);
        ++idx;
      }
    }
    return g;
  }

public:
  void set_kernel_args(const KernelArgs& args) override {
    if (args.grid_res > 0) grid_res_ = args.grid_res;
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("TPMKernel expects 3 parameter columns (cp, a0, b0), got %d",
                          (int)kernel_pars.cols.size());

             const int     n_comp  = static_cast<int>(comp_idx.size());
             const double* cp_col  = kernel_pars.cols[0];
             const double* a0_col  = kernel_pars.cols[1];
             const double* b0_col  = kernel_pars.cols[2];
             const double* cov_ptr = covariate.colptr(0);

             pred_mean_.resize(n_comp);
             pred_mode_.resize(n_comp);
             pred_logprecision_.resize(n_comp);

             const double cp_eps  = 1e-10;
             const TPMGrid grid   = build_grid(grid_res_);
             const int     nc     = grid.n_combi;
             const double  inv_nm1 = 1.0 / (nc - 1.0);

             std::vector<double> TPM_post(nc), TPM_pred(nc), TPM_update(nc);

             // initialise posterior with Beta prior from first trial
             {
               const int r0 = comp_idx[0];
               for (int k = 0; k < nc; ++k)
                 TPM_post[k] = dbeta_val(grid.p_XX[k], a0_col[r0], b0_col[r0])
                 * dbeta_val(grid.p_XY[k], a0_col[r0], b0_col[r0]);
               normalise_inplace(TPM_post);
             }

             for (int j = 0; j < n_comp; ++j) {
               const int    r       = comp_idx[j];
               const double cp      = cp_col[r];
               const double x       = cov_ptr[r];
               const bool   curr_na = is_nan(x);
               const bool   prev_na = (j == 0) || is_nan(cov_ptr[comp_idx[j - 1]]);
               const int    curr    = curr_na ? -1 : static_cast<int>(x);
               const int    prev    = (j == 0 || prev_na) ? -1
               : static_cast<int>(cov_ptr[comp_idx[j - 1]]);

               // degenerate: cp ≈ 1
               if ((1.0 - cp) < cp_eps) {
                 pred_mean_[j] = beta_mean(a0_col[r], b0_col[r]);
                 pred_mode_[j] = pred_mean_[j];
                 pred_logprecision_[j] = beta_log_precision(a0_col[r], b0_col[r]);
                 continue;
               }

               // degenerate: cp ≈ 0 — no volatility, read directly from posterior
               if (cp < cp_eps) {
                 pred_mean_[j] = prev_na
                 ? mean_discrete(grid.mean_p, TPM_post)
                   : (prev == 1 ? mean_discrete(grid.p_XX, TPM_post)
                        : mean_discrete(grid.p_XY, TPM_post));
                 pred_mode_[j] = pred_mean_[j];
                 pred_logprecision_[j] = prev_na
                 ? log_precision_discrete(grid.mean_p, TPM_post)
                   : (prev == 1 ? log_precision_discrete(grid.p_XX, TPM_post)
                        : log_precision_discrete(grid.p_XY, TPM_post));
               } else {
                 // full TPM update
                 const double sum_post = std::accumulate(
                   TPM_post.begin(), TPM_post.end(), 0.0);
                 const double mix_old  = 1.0 - cp;
                 const double mix_new  = cp;

                 for (int k = 0; k < nc; ++k)
                   TPM_pred[k] = mix_old * TPM_post[k]
                 + mix_new * (sum_post - TPM_post[k]) * inv_nm1;
                 normalise_inplace(TPM_pred);

                 pred_mean_[j] = prev_na
                 ? mean_discrete(grid.mean_p, TPM_pred)
                   : (prev == 1 ? mean_discrete(grid.p_XX, TPM_pred)
                        : mean_discrete(grid.p_XY, TPM_pred));
                 pred_mode_[j] = pred_mean_[j];
                 pred_logprecision_[j] = prev_na
                 ? log_precision_discrete(grid.mean_p, TPM_pred)
                   : (prev == 1 ? log_precision_discrete(grid.p_XX, TPM_pred)
                        : log_precision_discrete(grid.p_XY, TPM_pred));

                 if (curr_na || prev_na) {
                   TPM_post = TPM_pred;
                 } else {
                   const std::vector<double>* lp =
                     (prev == 0) ? (curr == 0 ? &grid.like_YY : &grid.like_XY)
                     : (curr == 0 ? &grid.like_YX : &grid.like_XX);
                   for (int k = 0; k < nc; ++k)
                     TPM_update[k] = mix_old * (*lp)[k] * TPM_post[k]
                   + mix_new * (*lp)[k]
                   * (sum_post - TPM_post[k]) * inv_nm1;
                   normalise_inplace(TPM_update);
                   std::swap(TPM_post, TPM_update);
                 }
               }
             }

             store_obs(cov_ptr, comp_idx);
             sync_out_to_mean();
             mark_run_complete();
           }
};


// ---- Type mapping + factory ----

KernelType to_kernel_type(const Rcpp::String& k);

std::unique_ptr<BaseKernel> make_kernel(KernelType kt,
                                        SEXP custom_fun = R_NilValue);


#endif // KERNELS_H

