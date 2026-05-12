#ifndef KERNELS_H
#define KERNELS_H

#include <unordered_map>
#include <memory>
#include <array>
#include <vector>     //
#include <Rcpp.h>    //
#include "nan_check.h"
#include "EMC2/userfun.hpp"

// View
struct KernelParsView {
  int n_rows;
  std::vector<const double*> cols;  // cols[k][row] = value for param k at trial row
};

// Struct for optional kernel arguments. Currently only "q-value resetting" is supported.
struct KernelArgs {
  const int* q_reset = nullptr;  // raw pointer into an IntegerVector; null = no reset
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
  RescorlaWagner
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
  // case KernelType::Delta2Kernel2:
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
  case KernelType::Custom: return{1, false};
  case KernelType::RescorlaWagner: return{-1, false};  // N columns allowed
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
                   const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& input,
           const std::vector<int>& comp_idx) override {

             const int n_comp   = static_cast<int>(comp_idx.size());
             const int n_pars   = static_cast<int>(kernel_pars.cols.size());
             const int n_inputs = input.ncol();

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
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);    // compressed output

             const double* cov_col = covariate.begin(); // column-major, col 0
             for (int j = 0; j < n_comp; ++j) {
               double x = cov_col[comp_idx[j]];
               out_[j] = is_nan(x) ? 0.0 : x;
             }

             mark_run_complete();
           }
};


struct LinDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col = covariate.begin(); // column-major, col 0
             for (int j = 0; j < n_comp; ++j) {
               double x = cov_col[comp_idx[j]];
               out_[j] = is_nan(x) ? 0.0 : -x;
             }

             mark_run_complete();
           }
};

struct ExpDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("ExpDecrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col    = covariate.begin();
             const double* lambda_col = kernel_pars.cols[0];
             for (int j = 0; j < n_comp; ++j) {
               int    r = comp_idx[j];
               double x = cov_col[r];
               out_[j] = is_nan(x) ? 0.0 : std::exp(-lambda_col[r] * x);
             }

             mark_run_complete();
           }
};

struct ExpIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("ExpIncrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col    = covariate.begin();
             const double* lambda_col = kernel_pars.cols[0];
             for (int j = 0; j < n_comp; ++j) {
               int    r = comp_idx[j];
               double x = cov_col[r];
               out_[j] = is_nan(x) ? 0.0 : 1.0 - std::exp(-lambda_col[r] * x);
             }


             mark_run_complete();
           }
};

struct PowDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("PowDecrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col   = covariate.begin();
             const double* alpha_col = kernel_pars.cols[0];
             for (int j = 0; j < n_comp; ++j) {
               int    r = comp_idx[j];
               double x = cov_col[r];
               out_[j] = is_nan(x) ? 0.0 : std::pow(1.0 + x, -alpha_col[r]);
             }

             mark_run_complete();
           }
};

struct PowIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("PowIncrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col   = covariate.begin();
             const double* alpha_col = kernel_pars.cols[0];
             for (int j = 0; j < n_comp; ++j) {
               int    r = comp_idx[j];
               double x = cov_col[r];
               out_[j] = is_nan(x) ? 0.0 : 1.0 - std::pow(1.0 + x, -alpha_col[r]);
             }

             mark_run_complete();
           }
};

struct Poly2Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("Poly2Kernel expects 2 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col = covariate.begin();
             const double* a1_col  = kernel_pars.cols[0];
             const double* a2_col  = kernel_pars.cols[1];
             for (int j = 0; j < n_comp; ++j) {
               int    r  = comp_idx[j];
               double x  = cov_col[r];
               if (!is_nan(x)) {
                 double x2 = x * x;
                 out_[j] = a1_col[r] * x + a2_col[r] * x2;
               }
             }

             mark_run_complete();
           }
};

struct Poly3Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("Poly3Kernel expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col = covariate.begin();
             const double* a1_col  = kernel_pars.cols[0];
             const double* a2_col  = kernel_pars.cols[1];
             const double* a3_col  = kernel_pars.cols[2];
             for (int j = 0; j < n_comp; ++j) {
               int    r  = comp_idx[j];
               double x  = cov_col[r];
               if (!is_nan(x)) {
                 double x2 = x * x;
                 out_[j] = a1_col[r] * x + a2_col[r] * x2 + a3_col[r] * x2 * x;
               }
             }

             mark_run_complete();
           }
};

struct Poly4Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 4) {
               Rcpp::stop("Poly4Kernel expects 4 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, 0);

             const double* cov_col = covariate.begin();
             const double* a1_col  = kernel_pars.cols[0];
             const double* a2_col  = kernel_pars.cols[1];
             const double* a3_col  = kernel_pars.cols[2];
             const double* a4_col  = kernel_pars.cols[3];
             for (int j = 0; j < n_comp; ++j) {
               int    r  = comp_idx[j];
               double x  = cov_col[r];
               if (!is_nan(x)) {
                 double x2 = x * x;
                 double x4 = x2 * x2;
                 out_[j] = a1_col[r] * x + a2_col[r] * x2 + a3_col[r] * x2 * x + a4_col[r] * x4;
               }
             }

             mark_run_complete();
           }
};


// Sequential kernels
struct SimpleDelta : DeltaKernel {
  SimpleDelta() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
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
             const double* cov_ptr   = covariate.begin();

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

struct Delta2LR : DeltaKernel {
  Delta2LR() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
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
// to-do - deprecate this, it doesn't work and is annoying to maintain...
struct Delta2Kernel : SequentialKernel {
  double qFast_ = NA_REAL;
  double qSlow_ = NA_REAL;
  double q_     = NA_REAL;
  const int* q_reset_ = nullptr;   // <-- ADD: null = no reset

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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("RescorlaWagnerKernel expects 2 parameter columns (q0, alpha), got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n_comp = static_cast<int>(comp_idx.size());
             n_covs_          = covariate.ncol();

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

// ---- Type mapping + factory ----

KernelType to_kernel_type(const Rcpp::String& k);

std::unique_ptr<BaseKernel> make_kernel(KernelType kt,
                                        SEXP custom_fun = R_NilValue);


#endif // KERNELS_H

