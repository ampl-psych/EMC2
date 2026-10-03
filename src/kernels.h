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
  int grid_res = 100;
  // const uint8_t* is_first_level_comp = nullptr;  // push mode only; null = filter mode
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
  SLinIncr,
  SLinDecr,
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
  case KernelType::SLinIncr:
  case KernelType::SLinDecr:
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
  int last_row_end = -1;  // -1 = never run

  // Remember expansion mapping (for 'at')
  // std::vector<int> expand_idx_;   // 1-based indices
  // bool has_expand_idx_ = false;

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
                   const std::vector<uint8_t>& at_mask,
                   const MatBool& nan_mask,
                   int row_start = 0,
                   int row_end   = -1) = 0;

  int rows_computed() const { return last_row_end; }

  virtual void reset() {
    last_row_end = -1;
    out_.clear();
    stream_buf_[0].clear();
    stream_buf_[1].clear();
    has_run_ = false;
  }

  // default no-op for non-sequential kernels
  virtual void rewind(int row_start) {}

  bool has_run() const { return has_run_; }

  bool has_run_for(int row_end_requested) const {
    return last_row_end == row_end_requested;
  }

  // const std::vector<double>& get_output() const { return out_; }

  // void set_expand_idx(const std::vector<int>& idx) {
  //   expand_idx_ = idx;
  //   has_expand_idx_ = !expand_idx_.empty();
  // }
  // const std::vector<int>& expand_idx() const { return expand_idx_; }
  // bool has_expand_idx() const { return has_expand_idx_; }

  // // Expand compressed out_ (length n_comp) into full length using expand_idx
  // void do_expand(const std::vector<int>& expand_idx) {
  //   const int n_full = static_cast<int>(expand_idx.size());
  //
  //   out_.resize(n_full);  // keeps compressed data in [0..n_comp-1]
  //
  //   for (int i = n_full - 1; i >= 0; --i) {
  //     int k = expand_idx[i] - 1;  // 1-based -> 0-based
  //     out_[i] = out_[k];
  //   }
  // }

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
  void mark_run_complete(int row_end) {
    has_run_     = true;
    last_row_end = row_end;
  }
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
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             const int n        = input.nrow;
             const int n_pars   = static_cast<int>(kernel_pars.cols.size());
             const int n_inputs = input.ncol;

             // collect active row indices from at_mask
             std::vector<int> active;
             active.reserve(n);
             for (int r = 0; r < n; ++r)
               if (at_mask[r]) active.push_back(r);

             const int n_active = static_cast<int>(active.size());
             out_.assign(n, 0.0);

             if (n_active == 0) { mark_run_complete(row_end < 0 ? n : row_end); return; }

             // build compressed parameter matrix: n_active x n_pars
             Rcpp::NumericMatrix pars_comp(n_active, n_pars);
             for (int p = 0; p < n_pars; ++p) {
               const double* col = kernel_pars.cols[p];
               for (int j = 0; j < n_active; ++j)
                 pars_comp(j, p) = col[active[j]];
             }

             // build compressed input matrix: n_active x n_inputs
             Rcpp::NumericMatrix input_comp(n_active, n_inputs);
             for (int j = 0; j < n_active; ++j) {
               const int r = active[j];
               for (int c = 0; c < n_inputs; ++c)
                 input_comp(j, c) = input(r, c);
             }

             // call user function
             userfun_t f = *fun_;
             Rcpp::NumericVector res = f(pars_comp, input_comp);

             if (res.size() != n_active)
               Rcpp::stop("CustomKernel: user function returned length %d, expected %d",
                          res.size(), n_active);

             // scatter back with carry-forward for non-active rows
             double last = 0.0;
             int    j    = 0;
             for (int r = 0; r < n; ++r) {
               if (at_mask[r]) last = res[j++];
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};

// For sequential kernels: currently same as BaseKernel -- just included to allow for other types (e.g. Bayesian ideal observer, autoregressive) in the future
struct SequentialKernel : BaseKernel {
  virtual ~SequentialKernel() {}
};

// All 1D delta kernels have scalar q and 1D pes_
struct DeltaKernel : SequentialKernel {
protected:
  double q_ = NA_REAL;
  double q_pending_ = NA_REAL;     // persists between incremental calls
  std::vector<double> pes_;        // PE per trial
  const int* q_reset_ = nullptr;   // null = no reset

public:
  virtual ~DeltaKernel() {}

  void set_kernel_args(const KernelArgs& args) override {
    q_reset_ = args.q_reset;
  }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 2);
  }

  KernelOutput get_output_stream(int code) const override {
    const int n = static_cast<int>(out_.size());
    if (code == 1) return KernelOutput{ out_.data(), n, 1 };
    if (code == 2) return KernelOutput{ pes_.data(), n, 1 };
    Rcpp::stop("DeltaKernel::get_output_stream: unsupported code %d (1=Q,2=PE)", code);
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qvalue";
    if (code == 2) return "PE";
    throw std::runtime_error("DeltaKernel::output_stream_name: unsupported code");
  }

  void reset() override {
    BaseKernel::reset();
    q_         = NA_REAL;
    q_pending_ = NA_REAL;
    pes_.clear();
  }

  void rewind(int row_start) override {
    q_         = out_[row_start];
    q_pending_ = out_[row_start];
  }
};


// ---- Individual kernels ----
// ---- Non-sequential kernels ----

struct LinIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) last = covariate(r, 0);
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};


struct LinDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) last = -covariate(r, 0);
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
            }
};

struct ExpDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("ExpDecrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* lambda_col = kernel_pars.cols[0];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) last = std::exp(-lambda_col[r] * covariate(r, 0));
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};

struct ExpIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("ExpIncrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* lambda_col = kernel_pars.cols[0];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) last = 1.0 - std::exp(-lambda_col[r] * covariate(r, 0));
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};

// Saturating linear kernels: slin_incr k = min(1, k_sat * c) rises linearly at
// rate k_sat and saturates at 1 once c >= 1/k_sat; slin_decr is its negative,
// k = -min(1, k_sat * c). Non-finite covariates (NA, NaN, Inf; e.g. SSD = Inf on
// go trials of a stop-signal design) give 0, so the trended parameter is left at
// its untrended value on those rows.
template <int SIGN>
struct SLinKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("SLinKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* k_col = kernel_pars.cols[0];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) {
                 const double x = covariate(r, 0);
                 if (is_finite(x)) {
                   const double v = k_col[r] * x;
                   last = SIGN * ((v < 1.0) ? v : 1.0);
                 } else {
                   last = 0.0;
                 }
               }
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};
using SLinIncrKernel = SLinKernel<1>;
using SLinDecrKernel = SLinKernel<-1>;

struct PowDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("PowDecrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* alpha_col = kernel_pars.cols[0];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) last = std::pow(1.0 + covariate(r, 0), -alpha_col[r]);
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};

struct PowIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 1) {
               Rcpp::stop("PowIncrKernel expects 1 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* alpha_col = kernel_pars.cols[0];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) last = 1.0 - std::pow(1.0 + covariate(r, 0), -alpha_col[r]);
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};

struct Poly2Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {
             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("Poly2Kernel expects 2 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) {
                 const double x = covariate(r, 0);
                 last = a1_col[r] * x + a2_col[r] * x * x;
               }
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};

struct Poly3Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("Poly3Kernel expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];
             const double* a3_col = kernel_pars.cols[2];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) {
                 const double x  = covariate(r, 0);
                 const double x2 = x * x;
                 last = a1_col[r] * x + a2_col[r] * x2 + a3_col[r] * x2 * x;
               }
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};

struct Poly4Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {
             if (kernel_pars.cols.size() != 4) {
               Rcpp::stop("Poly4Kernel expects 4 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n   = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             out_.resize(n);
             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];
             const double* a3_col = kernel_pars.cols[2];
             const double* a4_col = kernel_pars.cols[3];
             double last = (row_start > 0) ? out_[row_start - 1] : 0.0;
             for (int r = row_start; r < end; ++r) {
               if (at_mask[r]) {
                 const double x  = covariate(r, 0);
                 const double x2 = x * x;
                 last = a1_col[r] * x + a2_col[r] * x2
                 + a3_col[r] * x2 * x + a4_col[r] * x2 * x2;
               }
               out_[r] = last;
             }

             mark_run_complete(row_end < 0 ? n : row_end);
           }
};


// Sequential kernels
struct SimpleDelta : DeltaKernel {
  SimpleDelta() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {
             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("SimpleDelta expects 2 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }


             const int n = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             if (n <= 0) { out_.clear(); pes_.clear(); mark_run_complete(end); return; }

             const double*  q0_col    = kernel_pars.cols[0];
             const double*  alpha_col = kernel_pars.cols[1];
             const double*  cov_ptr   = covariate.colptr(0);
             const uint8_t* nm        = nan_mask.colptr(0);

             // initialise state only on first call
             if (row_start == 0) {
               q_         = q0_col[0];
               q_pending_ = q_;
               out_.resize(n);
               pes_.assign(n, NA_REAL);
             }
             // else: q_ and q_pending_ carry over from previous call

             for(int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // at_mask controls "commit" to pending Q-value
                 // ie., at the first level of `at`
                 q_ = q_pending_;
                 if (q_reset_ && q_reset_[r]) q_ = q0_col[r];
               }

               // Q-value is written on all rows
               out_[r] = q_;

               if(nm[r]) {
                 // not-nan mask controls whether q_pending_ needs updating
                 const double pe = cov_ptr[r] - q_;
                 pes_[r]         = pe;
                 q_pending_      = q_ + alpha_col[r] * pe;
               }
             }

             mark_run_complete(end);
           }
};

// Delta rule reparametrised to decouple the movement towards the outcome from the decay towards 0
struct DeltaDecoupled : DeltaKernel {
  DeltaDecoupled() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("DeltaDecoupled expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }
             const int n = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             if (n <= 0) { out_.clear(); pes_.clear(); mark_run_complete(end); return; }

             const double* q0_col     = kernel_pars.cols[0];
             const double* alpha_col  = kernel_pars.cols[1];
             const double* lambda_col = kernel_pars.cols[2];
             const double* cov_ptr    = covariate.colptr(0);
             const uint8_t* nm        = nan_mask.colptr(0);

             if (row_start == 0) {
               q_         = q0_col[0];
               q_pending_ = q_;
               out_.resize(n);
               pes_.assign(n, NA_REAL);
             }

             for(int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // at_mask controls commit
                 q_ = q_pending_;
                 if (q_reset_ && q_reset_[r]) q_ = q0_col[r];
               }

               // all rows are written
               out_[r] = q_;
               if(nm[r]) {
                 // not-nan mask controls update
                 const double x = cov_ptr[r];
                 pes_[r]        = x - q_;
                 q_pending_     = q_ + alpha_col[r] * x - lambda_col[r] * q_;
               }
             }
             mark_run_complete(end);
           }
};

struct Delta2LR : DeltaKernel {
  Delta2LR() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("Delta2LR expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;
             if (n <= 0) { out_.clear(); pes_.clear(); mark_run_complete(end); return; }

             const double* q0_col       = kernel_pars.cols[0];
             const double* alphaPos_col = kernel_pars.cols[1];
             const double* alphaNeg_col = kernel_pars.cols[2];
             const double*  cov_ptr   = covariate.colptr(0);
             const uint8_t* nm          = nan_mask.colptr(0);

             if (row_start == 0) {
               q_         = q0_col[0];
               q_pending_ = q_;
               out_.resize(n);
               pes_.assign(n, NA_REAL);
             }

             for(int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // at_mask controls commit (i.e., first level)
                 q_ = q_pending_;
                 if (q_reset_ && q_reset_[r]) q_ = q0_col[r];
               }
               // all rows write out
               out_[r] = q_;
               if(nm[r]) {
                 // not-nan mask controls update
                 const double pe    = cov_ptr[r] - q_;
                 const double alpha = (pe > 0.0) ? alphaPos_col[r] : alphaNeg_col[r];
                 pes_[r]            = pe;
                 q_pending_         = q_ + alpha * pe;
               }
             }

             mark_run_complete(end);
           }
};

// 2D PE kernel: separate from DeltaKernel
struct Delta2Kernel : SequentialKernel {
  double qFast_    = NA_REAL;
  double qSlow_    = NA_REAL;
  double q_        = NA_REAL;
  double qFast_pending_ = NA_REAL;
  double qSlow_pending_ = NA_REAL;
  double q_pending_     = NA_REAL;
  const int* q_reset_   = nullptr;

  // [compressed trial][0 = fast PE, 1 = slow PE]
  std::vector<double> pes_fast_;
  std::vector<double> pes_slow_;
  std::vector<double> q_fast_;
  std::vector<double> q_slow_;

  void reset() override {
    BaseKernel::reset();
    qFast_ = qSlow_ = q_ = NA_REAL;
    qFast_pending_ = qSlow_pending_ = q_pending_ = NA_REAL;
    q_fast_.clear(); q_slow_.clear();
    pes_fast_.clear(); pes_slow_.clear();
  }

  void rewind(int row_start) override {
    qFast_         = q_fast_[row_start];
    qSlow_         = q_slow_[row_start];
    q_             = out_[row_start];
    qFast_pending_ = qFast_;
    qSlow_pending_ = qSlow_;
    q_pending_     = q_;
  }

  void set_kernel_args(const KernelArgs& args) override {
    q_reset_ = args.q_reset;
  }


  Delta2Kernel() {}

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {
             if (kernel_pars.cols.size() != 4) {
               Rcpp::stop("Delta2Kernel expects 4 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n = covariate.nrow;
             const int end = (row_end < 0) ? n : row_end;

             if (n <= 0) {
               out_.clear(); q_fast_.clear(); q_slow_.clear();
               pes_fast_.clear(); pes_slow_.clear();
               mark_run_complete(end); return;
             }


             const double*  q0_col        = kernel_pars.cols[0];
             const double*  alphaFast_col = kernel_pars.cols[1];
             const double*  propSlow_col  = kernel_pars.cols[2];
             const double*  dSwitch_col   = kernel_pars.cols[3];
             const double* cov_ptr    = covariate.colptr(0);
             const uint8_t* nm            = nan_mask.colptr(0);

             if (row_start == 0) {
               qFast_ = qSlow_ = q_ = q0_col[0];
               qFast_pending_ = qSlow_pending_ = q_pending_ = q_;
               out_.resize(n);
               q_fast_.resize(n);
               q_slow_.resize(n);
               pes_fast_.assign(n, NA_REAL);
               pes_slow_.assign(n, NA_REAL);
             }

             for(int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // control commit
                 qFast_ = qFast_pending_;
                 qSlow_ = qSlow_pending_;
                 q_     = q_pending_;
                 if (q_reset_ && q_reset_[r])
                   qFast_ = qSlow_ = q_ = q0_col[r];
               }

               // always write out
               out_[r]    = q_;
               q_fast_[r] = qFast_;
               q_slow_[r] = qSlow_;

               if(nm[r]) {
                 // control update
                 const double x         = cov_ptr[r];
                 const double alphaFast = alphaFast_col[r];
                 const double alphaSlow = propSlow_col[r] * alphaFast;
                 const double dSwitch   = dSwitch_col[r];
                 const double peFast    = x - qFast_;
                 const double peSlow    = x - qSlow_;
                 pes_fast_[r]           = peFast;
                 pes_slow_[r]           = peSlow;
                 qFast_pending_          = qFast_ + alphaFast * peFast;
                 qSlow_pending_          = qSlow_ + alphaSlow * peSlow;
                 q_pending_              = (std::abs(qFast_pending_ - qSlow_pending_) > dSwitch) ? qFast_pending_ : qSlow_pending_;
               }
             }

             mark_run_complete(end);
           }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 5);
  }

  KernelOutput get_output_stream(int code) const override {
    const int n = static_cast<int>(out_.size());
    if (code == 1) return KernelOutput{ out_.data(), n, 1 };

    const std::vector<double>* src = nullptr;
    if      (code == 2) src = &q_fast_;
    else if (code == 3) src = &q_slow_;
    else if (code == 4) src = &pes_fast_;
    else if (code == 5) src = &pes_slow_;
    else Rcpp::stop("Delta2Kernel::get_output_stream: unsupported code %d", code);

    // output is already full-length — direct view, no copy needed
    return KernelOutput{ src->data(), n, 1 };
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
  std::vector<double> q_cur_;      // state per covariate
  std::vector<double> q_pending_;  // pending state per covariate
  std::vector<double> q_mat_;   // [n_comp * n_covs_]: Q-value per trial per covariate
  std::vector<double> pe_mat_;  // [n_comp * n_covs_]: compound PE for active covariates, NA otherwise

  const int* q_reset_ = nullptr;

public:
  void set_kernel_args(const KernelArgs& args) override {
    q_reset_ = args.q_reset;
    // if (args.is_first_level_comp != nullptr)
    //   Rcpp::stop("RescorlaWagnerKernel does not support at_mode = 'push'.");
  }

  void reset() override {
    BaseKernel::reset();
    q_mat_.clear();
    pe_mat_.clear();
    q_cur_.clear();
    q_pending_.clear();
    n_covs_ = 0;
  }

  void rewind(int row_start) override {
    for (int c = 0; c < n_covs_; ++c)
      q_cur_[c] = q_pending_[c] = q_mat_[row_start * n_covs_ + c];
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 2) {
               Rcpp::stop("RescorlaWagnerKernel expects 2 parameter columns (q0, alpha), got %d",
                          (int)kernel_pars.cols.size());
             }

             const int n = covariate.nrow;
             const int end   = (row_end < 0) ? n : row_end;
             const int ncov = covariate.ncol;

             if (n == 0 || ncov == 0) {
               q_mat_.clear(); pe_mat_.clear();
               mark_run_complete(end); return;
             }

             const double* q0_col    = kernel_pars.cols[0];
             const double* alpha_col = kernel_pars.cols[1];

             if (row_start == 0) {
               n_covs_ = ncov;
               q_mat_.assign(n * ncov, NA_REAL);
               pe_mat_.assign(n * ncov, NA_REAL);
               q_cur_.assign(ncov, q0_col[0]);
               q_pending_.assign(ncov, q0_col[0]);
             }


             for(int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // at_mask controls committing
                 q_cur_ = q_pending_;
                 if (q_reset_ && q_reset_[r]) {
                   const double q0_r = q0_col[r];
                   for(int c = 0; c < n_covs_; ++c) q_cur_[c] = q0_r;
                 }
               }

               // always write out
               for(int c = 0; c < n_covs_; ++c) q_mat_[r * n_covs_ + c] = q_cur_[c];

               // accumulate compound Q over active (non-NaN) covariates
               double reward    = NA_REAL;
               double q_active  = 0.0;
               bool   any_active = false;

               for(int c = 0; c < n_covs_; ++c) {
                 if(nan_mask(r, c)) {
                   // not-nan_mask controls update
                   reward    = covariate(r, c);
                   q_active += q_cur_[c];
                   any_active = true;
                 }
               }

               if(any_active) {
                 const double alpha       = alpha_col[r];
                 const double compound_pe = reward - q_active;
                 for (int c = 0; c < n_covs_; ++c) {
                   if(nan_mask(r, c)) {
                     pe_mat_[r * n_covs_ + c] = compound_pe;
                     q_pending_[c] = q_cur_[c] + alpha * compound_pe;
                   } else {
                     q_pending_[c] = q_cur_[c];
                   }
                 }
               }
             }

             mark_run_complete(end);
           }

  bool has_output_stream(int code) const override {
    return (code == 1 || code == 2);
  }

  // Stream 1: Q-matrix (n_rows x n_covs_), column-major
  // Stream 2: PE-matrix (n_rows x n_covs_), column-major
  KernelOutput get_output_stream(int code) const override {
    if (code != 1 && code != 2)
      Rcpp::stop("RescorlaWagnerKernel::get_output_stream: unsupported code %d", code);

    const std::vector<double>& src = (code == 1) ? q_mat_ : pe_mat_;
    const int n = static_cast<int>(src.size()) / n_covs_;
    std::vector<double>& buf = stream_buf_[code - 1];

    // transpose row-major [n x n_covs_] to column-major [n_covs_ x n] for R
    buf.resize(n * n_covs_);
    for (int c = 0; c < n_covs_; ++c)
      for (int r = 0; r < n; ++r)
        buf[c * n + r] = src[r * n_covs_ + c];

    return KernelOutput{ buf.data(), n, n_covs_ };
  }


  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qmatrix";
    if (code == 2) return "PEmatrix";
    throw std::runtime_error("RescorlaWagnerKernel::output_stream_name: unsupported code");
  }
};


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
  const int* belief_reset_ = nullptr;             //
  std::vector<double> n_hit_history_;
  std::vector<double> n_trial_history_;

  // incremental state — members, persist between run() calls
  double n_hit_         = 0.0;
  double n_trial_       = 0.0;
  double n_hit_pending_   = 0.0;
  double n_trial_pending_ = 0.0;
  double last_mean_ = 0.0;
  double last_mode_ = 0.0;
  double last_lp_   = 0.0;

  void store_obs(const double* cov_ptr, int n, const std::vector<uint8_t>& at_mask) {
    comp_obs_.resize(n);
    for (int r = 0; r < n; ++r)
      comp_obs_[r] = at_mask[r] ? cov_ptr[r] : std::numeric_limits<double>::quiet_NaN();
  }

  void ensure_surprise() const {
    if (surprise_computed_) return;
    const int n = static_cast<int>(pred_mean_.size());
    surprise_.resize(n, std::numeric_limits<double>::quiet_NaN());
    for (int r = 0; r < n; ++r)
      if (!is_nan(comp_obs_[r]))
        surprise_[r] = shannon_surprise(pred_mean_[r], comp_obs_[r]);
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
    n_hit_ = n_trial_ = n_hit_pending_ = n_trial_pending_ = 0.0;
    last_mean_ = last_mode_ = last_lp_ = 0.0;
    n_hit_history_.clear();
    n_trial_history_.clear();
  }

  void rewind(int row_start) override {
    n_hit_           = n_hit_history_[row_start];
    n_trial_         = n_trial_history_[row_start];
    n_hit_pending_   = n_hit_;
    n_trial_pending_ = n_trial_;
  }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 4);
  }

  void set_kernel_args(const KernelArgs& args) override {
    belief_reset_ = args.belief_reset;
  }

  KernelOutput get_output_stream(int code) const override {
    if (code == 3) ensure_surprise();
    const std::vector<double>* src = nullptr;
    if      (code == 1) src = &pred_mean_;
    else if (code == 2) src = &pred_mode_;
    else if (code == 3) src = &surprise_;
    else if (code == 4) src = &pred_logprecision_;
    else Rcpp::stop("DBMBaseKernel::get_output_stream: unsupported code %d", code);
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
  // keep out_ in sync with pred_mean_ for BaseKernel::get_output_stream (code=1)
  void sync_out_to_mean() { out_ = pred_mean_; }
};

// =============================================================================
// BetaBinomialKernel  —  basic (no memory constraint)
// Parameters: a0, b0
// =============================================================================

struct BetaBinomialKernel : DBMBaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 2)
               Rcpp::stop("BetaBinomialKernel expects 2 parameter columns (a0, b0), got %d",
                          (int)kernel_pars.cols.size());

             const int     n       = covariate.nrow;
             const int     end     = (row_end < 0) ? n : row_end;
             const double* a0_col  = kernel_pars.cols[0];
             const double* b0_col  = kernel_pars.cols[1];
             const double* cov_ptr = covariate.colptr(0);
             const uint8_t* nm     = nan_mask.colptr(0);

             if (row_start == 0) {
               n_hit_ = n_trial_ = n_hit_pending_ = n_trial_pending_ = 0.0;
               last_mean_ = last_mode_ = last_lp_ = 0.0;
               pred_mean_.resize(n);
               pred_mode_.resize(n);
               pred_logprecision_.resize(n);
               comp_obs_.resize(n);
               surprise_computed_ = false;
               n_hit_history_.assign(n, 0.0);
               n_trial_history_.assign(n, 0.0);
             }

             for (int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // at_mask controls commit
                 n_hit_   = n_hit_pending_;
                 n_trial_ = n_trial_pending_;
                 n_hit_history_[r]   = n_hit_;    // store after commit
                 n_trial_history_[r] = n_trial_;

                 if(belief_reset_ && belief_reset_[r]) {
                   n_hit_ = 0.0; n_trial_ = 0.0;
                 }
                 const double a_t = a0_col[r] + n_hit_;
                 const double b_t = b0_col[r] + (n_trial_ - n_hit_);
                 last_mean_ = pred_mean_[r]         = beta_mean(a_t, b_t);
                 last_mode_ = pred_mode_[r]         = beta_mode(a_t, b_t);
                 last_lp_   = pred_logprecision_[r] = beta_log_precision(a_t, b_t);
               } else {
                 pred_mean_[r]         = last_mean_;
                 pred_mode_[r]         = last_mode_;
                 pred_logprecision_[r] = last_lp_;
               }

               if(nm[r]) {
                 // not-nan mask controls update
                 n_hit_pending_   += cov_ptr[r];
                 n_trial_pending_ += 1.0;
                 }
             }

             store_obs(cov_ptr, n, at_mask);
             sync_out_to_mean();
             mark_run_complete(end);
           }
};

// =============================================================================
// BetaBinomialDecayKernel  —  exponential decay on accumulated counts
// Parameters: a0, b0, decay
// =============================================================================

struct BetaBinomialDecayKernel : DBMBaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("BetaBinomialDecayKernel expects 3 parameter columns "
                            "(a0, b0, decay), got %d",
                            (int)kernel_pars.cols.size());

             const int     n         = covariate.nrow;
             const int     end       = (row_end < 0) ? n : row_end;
             const double* a0_col    = kernel_pars.cols[0];
             const double* b0_col    = kernel_pars.cols[1];
             const double* decay_col = kernel_pars.cols[2];
             const double* cov_ptr   = covariate.colptr(0);
             const uint8_t* nm       = nan_mask.colptr(0);

             if (row_start == 0) {
               n_hit_ = n_trial_ = n_hit_pending_ = n_trial_pending_ = 0.0;
               last_mean_ = last_mode_ = last_lp_ = 0.0;
               pred_mean_.resize(n);
               pred_mode_.resize(n);
               pred_logprecision_.resize(n);
               comp_obs_.resize(n);
               surprise_computed_ = false;
               n_hit_history_.assign(n, 0.0);
               n_trial_history_.assign(n, 0.0);
             }

             for (int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // control commit
                 // apply one step of decay to pending counts (one tick per trial)
                 const double df = std::exp(-1.0 / decay_col[r]);
                 n_hit_pending_   = df * n_hit_pending_;
                 n_trial_pending_ = df * n_trial_pending_;

                 // then commit
                 n_hit_   = n_hit_pending_;
                 n_trial_ = n_trial_pending_;
                 n_hit_history_[r]   = n_hit_;    // store after commit
                 n_trial_history_[r] = n_trial_;

                 if (belief_reset_ && belief_reset_[r]) {
                   n_hit_ = n_hit_pending_ = 0.0;
                   n_trial_ = n_trial_pending_ = 0.0;
                 }
                 const double a_t = a0_col[r] + n_hit_;
                 const double b_t = b0_col[r] + (n_trial_ - n_hit_);
                 last_mean_ = pred_mean_[r]         = beta_mean(a_t, b_t);
                 last_mode_ = pred_mode_[r]         = beta_mode(a_t, b_t);
                 last_lp_   = pred_logprecision_[r] = beta_log_precision(a_t, b_t);
               } else {
                 pred_mean_[r]         = last_mean_;
                 pred_mode_[r]         = last_mode_;
                 pred_logprecision_[r] = last_lp_;
               }

               if(nm[r]) {
                 // trigger update
                 n_hit_pending_   += cov_ptr[r];
                 n_trial_pending_ += 1.0;
               }
             }

             store_obs(cov_ptr, n, at_mask);
             sync_out_to_mean();
             mark_run_complete(end);
           }
};

// =============================================================================
// BetaBinomialWindowKernel  —  fixed sliding window
// Parameters: a0, b0, window
// =============================================================================

struct BetaBinomialWindowKernel : DBMBaseKernel {
private:
  struct Event { double obs; int idx; };
  std::deque<Event> buf_;
  std::deque<Event> buf_pending_;
  std::vector<std::deque<Event>> buf_history_;

public:
  void reset() override {
    DBMBaseKernel::reset();
    buf_.clear();
    buf_pending_.clear();
    buf_history_.clear();
  }

  void rewind(int row_start) override {
    DBMBaseKernel::rewind(row_start);
    buf_         = buf_history_[row_start];
    buf_pending_ = buf_;
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("BetaBinomialWindowKernel expects 3 parameter columns "
                            "(a0, b0, window), got %d",
                            (int)kernel_pars.cols.size());

             const int     n          = covariate.nrow;
             const int     end        = (row_end < 0) ? n : row_end;
             const double* a0_col     = kernel_pars.cols[0];
             const double* b0_col     = kernel_pars.cols[1];
             const double* window_col = kernel_pars.cols[2];
             const double* cov_ptr    = covariate.colptr(0);
             const uint8_t* nm        = nan_mask.colptr(0);

             if (row_start == 0) {
               n_hit_ = n_trial_ = n_hit_pending_ = n_trial_pending_ = 0.0;
               last_mean_ = last_mode_ = last_lp_ = 0.0;
               buf_.clear(); buf_pending_.clear();
               pred_mean_.resize(n);
               pred_mode_.resize(n);
               pred_logprecision_.resize(n);
               comp_obs_.resize(n);
               surprise_computed_ = false;
               n_hit_history_.assign(n, 0.0);
               n_trial_history_.assign(n, 0.0);
               buf_history_.assign(n, std::deque<Event>{});
             }

             for (int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // commit
                 n_hit_   = n_hit_pending_;
                 n_trial_ = n_trial_pending_;
                 buf_     = buf_pending_;

                 buf_history_[r]         = buf_;
                 n_hit_history_[r]       = n_hit_;
                 n_trial_history_[r]     = n_trial_;

                 if (belief_reset_ && belief_reset_[r]) {
                   n_hit_ = 0.0; n_trial_ = 0.0; buf_.clear();
                 }

                 const int w = static_cast<int>(window_col[r]);
                 while (!buf_.empty() && (r - buf_.front().idx) > w) {
                   n_hit_   -= buf_.front().obs;
                   n_trial_ -= 1.0;
                   buf_.pop_front();
                 }

                 // snapshot pending from newly committed state
                 buf_pending_     = buf_;
                 n_hit_pending_   = n_hit_;
                 n_trial_pending_ = n_trial_;

                 const double a_t = a0_col[r] + n_hit_;
                 const double b_t = b0_col[r] + (n_trial_ - n_hit_);
                 last_mean_ = pred_mean_[r]         = beta_mean(a_t, b_t);
                 last_mode_ = pred_mode_[r]         = beta_mode(a_t, b_t);
                 last_lp_   = pred_logprecision_[r] = beta_log_precision(a_t, b_t);
               } else {
                 pred_mean_[r]         = last_mean_;
                 pred_mode_[r]         = last_mode_;
                 pred_logprecision_[r] = last_lp_;
               }

               if(nm[r]) {
                 buf_pending_.push_back({cov_ptr[r], r});
                 n_hit_pending_   += cov_ptr[r];
                 n_trial_pending_ += 1.0;
               }
             }

             store_obs(cov_ptr, n, at_mask);
             sync_out_to_mean();
             mark_run_complete(end);
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
  std::vector<double> DBM_post_;
  std::vector<double> DBM_post_pending_;
  bool first_active_ = true;
  std::vector<std::vector<double>> post_history_;

public:
  void reset() override {
    DBMBaseKernel::reset();
    DBM_post_.clear();
    DBM_post_pending_.clear();
    first_active_ = true;
    post_history_.clear();
  }

  void rewind(int row_start) override {
    DBMBaseKernel::rewind(row_start);
    DBM_post_         = post_history_[row_start];
    DBM_post_pending_ = DBM_post_;
    first_active_     = false;  // we've been here before
  }

  void set_kernel_args(const KernelArgs& args) override {
    DBMBaseKernel::set_kernel_args(args);
    if (args.grid_res > 0) grid_res_ = args.grid_res;
    // if (args.is_first_level_comp != nullptr)
    //   Rcpp::stop("DBMKernel does not support at_mode = 'push'.");
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("DBMKernel expects 3 parameter columns (cp, mu0, s0), got %d",
                          (int)kernel_pars.cols.size());

             const int     n       = covariate.nrow;
             const int     end     = (row_end < 0) ? n : row_end;
             const double* cp_col  = kernel_pars.cols[0];
             const double* mu0_col = kernel_pars.cols[1];
             const double* s0_col  = kernel_pars.cols[2];
             const double* cov_ptr = covariate.colptr(0);
             const uint8_t* nm     = nan_mask.colptr(0);

             const int gs = grid_res_ + 1;
             std::vector<double> prob_grid(gs), x_like(gs), y_like(gs);
             for (int i = 0; i < gs; ++i) {
               prob_grid[i] = static_cast<double>(i) / (gs - 1);
               x_like[i]   = prob_grid[i];
               y_like[i]   = 1.0 - prob_grid[i];
             }

             if (row_start == 0) {
               DBM_post_.assign(gs, 0.0);
               DBM_post_pending_.assign(gs, 0.0);
               first_active_ = true;
               last_mean_ = last_mode_ = last_lp_ = 0.0;
               pred_mean_.resize(n);
               pred_mode_.resize(n);
               pred_logprecision_.resize(n);
               comp_obs_.resize(n);
               surprise_computed_ = false;
               post_history_.assign(n, std::vector<double>{});
             }
             for (int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 DBM_post_ = DBM_post_pending_;
                 post_history_[r] = DBM_post_;  // store committed posterior

                 const double cp  = cp_col[r];
                 const double mu0 = mu0_col[r];
                 const double s0  = s0_col[r];
                 const double a   = mu0 * s0;
                 const double b   = (1.0 - mu0) * s0;
                 const bool reset = first_active_ || (belief_reset_ && belief_reset_[r]);

                 std::vector<double> DBM_prior(gs), DBM_pred(gs);
                 for (int i = 0; i < gs; ++i) DBM_prior[i] = dbeta_val(prob_grid[i], a, b);
                 normalise_inplace(DBM_prior);

                 if (reset) {
                   DBM_pred = DBM_prior;
                 } else {
                   for (int i = 0; i < gs; ++i) DBM_pred[i] = (1.0 - cp) * DBM_post_[i] + cp * DBM_prior[i];
                   normalise_inplace(DBM_pred);
                 }

                 last_mean_ = pred_mean_[r]         = mean_discrete(prob_grid, DBM_pred);
                 last_mode_ = pred_mode_[r]         = mode_discrete(prob_grid, DBM_pred);
                 last_lp_   = pred_logprecision_[r] = log_precision_discrete(prob_grid, DBM_pred);

                 // compute pending posterior
                 if (!nm[r]) {
                   DBM_post_pending_ = DBM_pred;
                 } else {
                   const double x = cov_ptr[r];
                   const std::vector<double>& like = (x == 1.0) ? x_like : y_like;
                   for (int i = 0; i < gs; ++i) DBM_post_pending_[i] = DBM_pred[i] * like[i];
                   normalise_inplace(DBM_post_pending_);
                 }

                 first_active_ = false;
               } else {
                 pred_mean_[r]         = last_mean_;
                 pred_mode_[r]         = last_mode_;
                 pred_logprecision_[r] = last_lp_;
               }
             }
             store_obs(cov_ptr, n, at_mask);
             sync_out_to_mean();
             mark_run_complete(end);
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
  std::vector<double> TPM_post_;
  std::vector<double> TPM_post_pending_;
  bool first_active_  = true;
  int  prev_active_r_ = -1;
  std::vector<std::vector<double>> post_history_;
  std::vector<int>                 prev_active_r_history_;

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
  void reset() override {
    DBMBaseKernel::reset();
    TPM_post_.clear();
    TPM_post_pending_.clear();
    first_active_  = true;
    prev_active_r_ = -1;
    post_history_.clear();
    prev_active_r_history_.clear();
  }

  void rewind(int row_start) override {
    DBMBaseKernel::rewind(row_start);
    TPM_post_         = post_history_[row_start];
    TPM_post_pending_ = TPM_post_;
    prev_active_r_    = prev_active_r_history_[row_start];
    first_active_     = false;
  }

  void set_kernel_args(const KernelArgs& args) override {
    DBMBaseKernel::set_kernel_args(args);
    if (args.grid_res > 0) grid_res_ = args.grid_res;
    // if (args.is_first_level_comp != nullptr)
    //   Rcpp::stop("TPMKernel does not support at_mode = 'push'.");
  }

  void run(const KernelParsView& kernel_pars,
           const Mat& covariate,
           const std::vector<uint8_t>& at_mask,
           const MatBool& nan_mask,
           int row_start = 0,
           int row_end   = -1) override {

             if (kernel_pars.cols.size() != 3)
               Rcpp::stop("TPMKernel expects 3 parameter columns (cp, a0, b0), got %d",
                          (int)kernel_pars.cols.size());

             const int     n       = covariate.nrow;
             const int     end     = (row_end < 0) ? n : row_end;
             const double* cp_col  = kernel_pars.cols[0];
             const double* a0_col  = kernel_pars.cols[1];
             const double* b0_col  = kernel_pars.cols[2];
             const double* cov_ptr = covariate.colptr(0);
             const uint8_t* nm     = nan_mask.colptr(0);

             const TPMGrid  grid    = build_grid(grid_res_);
             const int      nc      = grid.n_combi;
             const double   inv_nm1 = 1.0 / (nc - 1.0);

             if (row_start == 0) {
               TPM_post_.assign(nc, 0.0);
               TPM_post_pending_.assign(nc, 0.0);
               first_active_  = true;
               prev_active_r_ = -1;
               last_mean_ = last_mode_ = last_lp_ = 0.0;
               pred_mean_.resize(n);
               pred_mode_.resize(n);
               pred_logprecision_.resize(n);
               comp_obs_.resize(n);
               surprise_computed_ = false;
               post_history_.assign(n, std::vector<double>{});
               prev_active_r_history_.assign(n, -1);
             }
             std::vector<double> TPM_pred(nc), TPM_update(nc);

             for(int r = row_start; r < end; ++r) {
               if (at_mask[r]) {
                 TPM_post_ = TPM_post_pending_;
                 post_history_[r]          = TPM_post_;
                 prev_active_r_history_[r] = prev_active_r_;  // store *before* updating prev_active_r_

                 const double cp      = cp_col[r];
                 const double x       = cov_ptr[r];
                 const bool   reset   = first_active_ || (belief_reset_ && belief_reset_[r]);
                 const bool   prev_na = (prev_active_r_ < 0) || !nm[prev_active_r_];
                 const int    prev    = prev_na ? -1 : static_cast<int>(cov_ptr[prev_active_r_]);
                 const bool   curr_na = !nm[r];
                 const int    curr    = curr_na ? -1 : static_cast<int>(x);

                 const double sum_post = std::accumulate(TPM_post_.begin(), TPM_post_.end(), 0.0);

                 if(reset) {
                   for(int k = 0; k < nc; ++k) TPM_pred[k] = dbeta_val(grid.p_XX[k], a0_col[r], b0_col[r]) * dbeta_val(grid.p_XY[k], a0_col[r], b0_col[r]);
                 } else {
                   for(int k = 0; k < nc; ++k) TPM_pred[k] = (1.0 - cp) * TPM_post_[k] + cp * (sum_post - TPM_post_[k]) * inv_nm1;
                 }
                 normalise_inplace(TPM_pred);

                 if(prev_na || reset) {
                   last_mean_ = pred_mean_[r]         = mean_discrete(grid.mean_p, TPM_pred);
                   last_mode_ = pred_mode_[r]         = mode_discrete(grid.mean_p, TPM_pred);
                   last_lp_   = pred_logprecision_[r] = log_precision_discrete(grid.mean_p, TPM_pred);
                 } else {
                   last_mean_ = pred_mean_[r]         = prev == 1 ? mean_discrete(grid.p_XX, TPM_pred) : mean_discrete(grid.p_XY, TPM_pred);
                   last_mode_ = pred_mode_[r]         = prev == 1 ? mode_discrete(grid.p_XX, TPM_pred) : mode_discrete(grid.p_XY, TPM_pred);
                   last_lp_   = pred_logprecision_[r] = prev == 1 ? log_precision_discrete(grid.p_XX, TPM_pred) : log_precision_discrete(grid.p_XY, TPM_pred);
                 }

                 // compute pending posterior
                 if (curr_na || prev_na || reset) {
                   TPM_post_pending_ = TPM_pred;
                 } else {
                   const std::vector<double>* lp = (prev == 0) ? (curr == 0 ? &grid.like_YY : &grid.like_XY) : (curr == 0 ? &grid.like_YX : &grid.like_XX);
                   for (int k = 0; k < nc; ++k) TPM_update[k] = (1.0 - cp) * (*lp)[k] * TPM_post_[k] + cp * (*lp)[k] * (sum_post - TPM_post_[k]) * inv_nm1;
                   normalise_inplace(TPM_update);
                   std::swap(TPM_post_pending_, TPM_update);
                 }

                 first_active_  = false;
                 prev_active_r_ = r;
               } else {
                 pred_mean_[r]         = last_mean_;
                 pred_mode_[r]         = last_mode_;
                 pred_logprecision_[r] = last_lp_;
               }
             }

             store_obs(cov_ptr, n, at_mask);
             sync_out_to_mean();
             mark_run_complete(end);
           }
};


// ---- Type mapping + factory ----

KernelType to_kernel_type(const Rcpp::String& k);

std::unique_ptr<BaseKernel> make_kernel(KernelType kt,
                                        SEXP custom_fun = R_NilValue);


#endif // KERNELS_H

