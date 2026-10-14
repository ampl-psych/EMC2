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
    return has_run_ && last_row_end == row_end_requested;
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
    has_run_ = false;
    last_row_end = row_start;
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
  std::vector<double> surprise_;
  std::vector<double> pred_logprecision_;
  const int* belief_reset_ = nullptr;
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

  void init_streams(int n) {
    last_mean_ = last_mode_ = last_lp_ = 0.0;
    pred_mean_.resize(n);
    pred_mode_.resize(n);
    pred_logprecision_.resize(n);
    surprise_.assign(n, NA_REAL);
  }

  // active row: new prediction; surprise defaults to NA until an observation arrives
  void set_active(int r, double mean, double mode, double lp) {
    last_mean_ = pred_mean_[r]         = mean;
    last_mode_ = pred_mode_[r]         = mode;
    last_lp_   = pred_logprecision_[r] = lp;
    surprise_[r] = NA_REAL;
  }

  // inactive row: carry forward, no surprise
  void set_inactive(int r) {
    pred_mean_[r]         = last_mean_;
    pred_mode_[r]         = last_mode_;
    pred_logprecision_[r] = last_lp_;
    surprise_[r]          = NA_REAL;
  }

  // call only on active rows with a valid observation
  void set_surprise(int r, double obs) {
    surprise_[r] = shannon_surprise(last_mean_, obs);
  }

public:
  void reset() override {
    BaseKernel::reset();
    pred_mean_.clear();
    pred_mode_.clear();
    surprise_.clear();
    pred_logprecision_.clear();
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

    has_run_ = false;
    last_row_end = row_start;
  }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 4);
  }

  void set_kernel_args(const KernelArgs& args) override {
    belief_reset_ = args.belief_reset;
  }

  KernelOutput get_output_stream(int code) const override {
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
               init_streams(n);
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
                   n_hit_pending_ = 0.0; n_trial_pending_ = 0.0;
                 }
                 const double a_t = a0_col[r] + n_hit_;
                 const double b_t = b0_col[r] + (n_trial_ - n_hit_);
                 double mean, mode, log_precision;
                 beta_summary(a_t, b_t, mean, mode, log_precision); // inline calculated
                 set_active(r, mean, mode, log_precision);
               } else {
                 set_inactive(r);
               }

               if(nm[r]) {
                 // not-nan mask controls update
                 set_surprise(r, cov_ptr[r]);
                 n_hit_pending_   += cov_ptr[r];
                 n_trial_pending_ += 1.0;
               }
             }
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
               init_streams(n);
               n_hit_history_.assign(n, 0.0);
               n_trial_history_.assign(n, 0.0);
             }

             for (int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 n_hit_history_[r]   = n_hit_pending_;    // store after commit
                 n_trial_history_[r] = n_trial_pending_;

                 // control commit
                 // apply one step of decay to pending counts (one tick per trial)
                 const double df = std::exp(-1.0 / decay_col[r]);
                 n_hit_pending_   = df * n_hit_pending_;
                 n_trial_pending_ = df * n_trial_pending_;

                 // then commit
                 n_hit_   = n_hit_pending_;
                 n_trial_ = n_trial_pending_;

                 if (belief_reset_ && belief_reset_[r]) {
                   n_hit_ = n_hit_pending_ = 0.0;
                   n_trial_ = n_trial_pending_ = 0.0;
                 }
                 const double a_t = a0_col[r] + n_hit_;
                 const double b_t = b0_col[r] + (n_trial_ - n_hit_);
                 double mean, mode, log_precision;
                 beta_summary(a_t, b_t, mean, mode, log_precision); // inline calculated
                 set_active(r, mean, mode, log_precision);
               } else {
                set_inactive(r);
              }

               if(nm[r]) {
                 // trigger update
                 set_surprise(r, cov_ptr[r]);
                 n_hit_pending_   += cov_ptr[r];
                 n_trial_pending_ += 1.0;
               }
             }
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
               buf_.clear(); buf_pending_.clear();
               init_streams(n);
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
                 double mean, mode, log_precision;
                 beta_summary(a_t, b_t, mean, mode, log_precision); // inline calculated
                 set_active(r, mean, mode, log_precision);
               } else {
                 set_inactive(r);
               }

               if(nm[r]) {
                 buf_pending_.push_back({cov_ptr[r], r});
                 set_surprise(r, cov_ptr[r]);
                 n_hit_pending_   += cov_ptr[r];
                 n_trial_pending_ += 1.0;
               }
             }

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
  int grid_built_for_ = -1;
  std::vector<double> DBM_post_;
  std::vector<double> DBM_post_pending_;
  std::vector<std::vector<double>> post_history_;
  bool first_active_ = true;  // have we started yet?
  int first_trial_row_ = -1;


  std::vector<double> prob_grid_, x_like_, y_like_;
  // Scratch vectors to write in posterior probabilities across the grid
  // Resized in ensure_grid();
  // DBM_prior_ : cache of the normalised Beta(a, b) prior, keyed on (a, b); see ensure_prior()
  std::vector<double> DBM_prior_, DBM_pred_;

  std::vector<double> log_x_, log_1mx_;                 // grid-only log terms
  double cached_a_ = std::numeric_limits<double>::quiet_NaN();
  double cached_b_ = std::numeric_limits<double>::quiet_NaN();
  bool   prior_valid_ = false;

  // Build a grid once, and calculate the dbeta_values once per grid.
  // log p0(\gamma_i) = (a-1)*log(\gamma_i) + (b-1)*log(1-\gamma_i) - [ln\Gamma(a) + ln\Gamma(b) - ln\Gamma(a+b)]
  // log(\gamma_i) and log(1-\gamma_i) depend only on i - so only on the grid, not on a or b -- precompute them once rather than each row.
  // [ln\Gamma(a) + ln\Gamma(b) - ln\Gamma(a+b)]  is constant across i, so it cancels in the normalisation. It can be removed entirely.
  void ensure_grid() {
    if (grid_built_for_ == grid_res_) return;
    const int gs = grid_res_ + 1;
    prob_grid_.resize(gs); x_like_.resize(gs); y_like_.resize(gs);
    log_x_.assign(gs, 0.0); log_1mx_.assign(gs, 0.0);
    DBM_prior_.resize(gs); DBM_pred_.resize(gs);
    for (int i = 0; i < gs; ++i) {
      prob_grid_[i] = static_cast<double>(i) / (gs - 1);
      x_like_[i]    = prob_grid_[i];
      y_like_[i]    = 1.0 - prob_grid_[i];
    }
    for (int i = 1; i < gs - 1; ++i) {                  // endpoints unused
      log_x_[i]   = std::log(prob_grid_[i]);            // log(\gamma_i)
      log_1mx_[i] = std::log1p(-prob_grid_[i]);         // log(1-\gamma_i)
    }
    prior_valid_    = false;                            // grid changed -> cached prior is stale
    grid_built_for_ = grid_res_;
  }

  // DBM_prior_ is now a cache keyed on (a, b), not scratch space
  void ensure_prior(double a, double b) {
    if (prior_valid_ && a == cached_a_ && b == cached_b_) return;
    beta_prior_grid(log_x_, log_1mx_, a, b, DBM_prior_);  // <- calculates the Beta prior across the grid with a and b. See kernels_math.h
    cached_a_ = a; cached_b_ = b; prior_valid_ = true;    // keep this in a cache and only recalculate when the prior changes.
  }


public:
  void reset() override {
    DBMBaseKernel::reset();
    DBM_post_.clear();
    DBM_post_pending_.clear();
    first_active_ = true;
    first_trial_row_ = -1;
    post_history_.clear();
  }

  void rewind(int row_start) override {
    if(post_history_[row_start].empty()) Rcpp::stop("DBMKernel::rewind: row %d is not a trial start", row_start);
    DBM_post_pending_ = post_history_[row_start];
    first_active_ = (first_trial_row_ < 0) || (row_start <= first_trial_row_);

    has_run_ = false;
    last_row_end = row_start;
  }

  void set_kernel_args(const KernelArgs& args) override {
    DBMBaseKernel::set_kernel_args(args);
    if (args.grid_res > 0) grid_res_ = args.grid_res;
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

             // Grid is member of the kernel, no need to re-create it each run() call
             ensure_grid();
             const int gs = grid_res_ + 1;
             const std::vector<double>& prob_grid = prob_grid_;
             const std::vector<double>& x_like    = x_like_;
             const std::vector<double>& y_like    = y_like_;
             std::vector<double>&       DBM_prior = DBM_prior_;
             std::vector<double>&       DBM_pred  = DBM_pred_;

             if (row_start == 0) {
               DBM_post_pending_.assign(gs, 0.0);
               first_active_ = true;
               first_trial_row_ = -1;
               init_streams(n);
               post_history_.assign(n, std::vector<double>{});
             }

             for (int r = row_start; r < end; ++r) {
               if(at_mask[r]) {
                 // at mask controls committing to stored value
                 // Last trial's pending posterior is committed:
                 // DBM_post_[i]  =  p(\gamma_{t-1} = \gamma_i | x_{1:t-1})
                 post_history_[r] = DBM_post_pending_;  // store committed posterior, to be able to rewind

                 // Beta prior: p0 ~ Beta(a,b)
                 // a = \mu0*s0, b = (1-\mu0)·s0
                 // DBM_prior[i] = p0(\gamma_i), normalised so that \Sum_i p0_i = 1
                 const double cp  = cp_col[r];
                 const double mu0 = mu0_col[r];
                 const double s0  = s0_col[r];
                 const double a   = mu0 * s0;
                 const double b   = (1.0 - mu0) * s0;

                 // Prior probability across grid values:
                 // Prior p0(\gamma_i) = Beta(a, b) on the grid, normalised (cached on (a, b))
                 ensure_prior(a, b);

                 const bool reset = first_active_ || (belief_reset_ && belief_reset_[r]);
                 if(reset) {
                   // No memory (first trial or belief reset): pred = p0, i.e. c = 1
                   DBM_pred = DBM_prior;
                 } else {
                   //  --- TRANSITION -----------------------------------------------
                   // The probability of gamma is a mixture: Either it's drawn from the prior (with probability c)
                   // or it's the same as last (with probability 1-c):
                   // p(\gamma_t | x_{1:t-1}) = (1-c)*p(\gamma_{t-1} | x_{1:t-1}) + c*p0(\gamma_t)
                   // DBM_pred[i]             = (1-c)*DBM_post_[i]                + c*DBM_prior[i]
                   const double keep = 1.0 - cp;
#pragma omp simd
                   for(int i = 0; i < gs; ++i) DBM_pred[i] = keep * DBM_post_pending_[i] + cp * DBM_prior[i];
                   normalise_inplace(DBM_pred);  // ensure integration to one
                 }

                 // --- PREDICTION -----------------------------------------------
                 // Summaries of p(\gamma_t | x_{1:t-1}):
                 //   mean = \Sum_i \gamma_i*w_i  = P(x_t = 1 | x_{1:t-1})   (the prediction)
                 //   mode = \gamma_i at the largest w_i
                 //   log-precision = -log Var(\gamma_t),  Var = \Sum_i (\gamma_i - mean)**2 * w_i
                 double mean, mode, log_precision;
                 discrete_summary(prob_grid, DBM_pred, mean, mode, log_precision);  // calculated in-place; shares arithmic
                 set_active(r, mean, mode, log_precision);

                 // Carry forward current predictions until next observation nm
                 DBM_post_pending_ = DBM_pred;

                 // We've started
                 if(first_active_) first_trial_row_ = r;
                 first_active_ = false;
               } else {
                 set_inactive(r);  // carry forward current beliefs
               }

               // ========================================================================
               // OBSERVATION (nm): first non-NaN value of the trial, on any row
               // ========================================================================
               if(nm[r]) {
                 // --- SURPRISE ---------------------------------------------------------
                 //     S_t = -log2 P(x_t | x_{1:t-1})
                 //         = -log2(\hat{p}_t)      if x = 1
                 //         = -log2(1 - \hat{p}_t)  if x = 0, with \hat{p}_t = last_mean_
                 const double x = cov_ptr[r];
                 set_surprise(r, x);

                 // --- BAYESIAN UPDATE ------------------------------------------
                 //     p(\gamma_t | x_{1:t}) \propto \gamma^x (1-\gamma)^(1-x) * p(\gamma_t | x_{1:t-1})
                 //     like[i] = \gamma_i      if x = 1   (x_like)
                 //             = 1 - \gamma_i  if x = 0   (y_like)
                 //     DBM_post_pending_ currently holds DBM_pred, so multiply in place
                 //     and renormalise (the division by P(x_t | x_{1:t-1})).
                 const std::vector<double>& like = (x == 1.0) ? x_like : y_like;
#pragma omp simd
                 for (int i = 0; i < gs; ++i) DBM_post_pending_[i] *= like[i];
                 normalise_inplace(DBM_post_pending_);
               }
             }
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

  int grid_built_for_ = -1;
  TPMGrid grid_;
  std::vector<double> TPM_pred_;
  std::vector<double> TPM_prior_;
  std::vector<double> beta_1d_;

  double cached_a_ = std::numeric_limits<double>::quiet_NaN();
  double cached_b_ = std::numeric_limits<double>::quiet_NaN();
  bool prior_valid_ = false;

  // Create grid
  void ensure_grid() {
    if (grid_built_for_ == grid_res_) return;

    grid_ = build_grid(grid_res_);
    TPM_pred_.resize(grid_.n_combi);
    TPM_prior_.resize(grid_.n_combi);
    beta_1d_.resize(grid_.resol);

    prior_valid_ = false;
    grid_built_for_ = grid_res_;
  }

  // Create prior -- calculate on grid where possible, and keep in cache.
  void ensure_prior(double a, double b) {
    if (prior_valid_ && a == cached_a_ && b == cached_b_) return;

    const int resol = grid_.resol;

    // The first row of p_XX contains the complete 1D grid.
    for (int i = 0; i < resol; ++i) {
      beta_1d_[i] = dbeta_val(grid_.p_XX[i], a, b);
    }

    for (int i0 = 0; i0 < resol; ++i0) {
      const double beta_xy = beta_1d_[i0];
      const int offset = i0 * resol;

#pragma omp simd
      for (int i1 = 0; i1 < resol; ++i1) {
        TPM_prior_[offset + i1] = beta_1d_[i1] * beta_xy;
      }
    }

    normalise_inplace(TPM_prior_);

    cached_a_ = a;
    cached_b_ = b;
    prior_valid_ = true;
  }

public:
  void reset() override {
    DBMBaseKernel::reset();
    TPM_post_pending_.clear();
    first_active_  = true;
    prev_active_r_ = -1;
    post_history_.clear();
    prev_active_r_history_.clear();
  }

  void rewind(int row_start) override {
    if(row_start < 0 || static_cast<std::size_t>(row_start) >= post_history_.size() || post_history_[row_start].empty()) {
      Rcpp::stop( "TPMKernel::rewind: row %d is not a recorded trial start", row_start);
    }

    TPM_post_pending_ = post_history_[row_start];
    prev_active_r_ = prev_active_r_history_[row_start];
    first_active_ = (prev_active_r_ < 0);

    has_run_ = false;
    last_row_end = row_start;

  }

  void set_kernel_args(const KernelArgs& args) override {
    DBMBaseKernel::set_kernel_args(args);
    if (args.grid_res > 0) grid_res_ = args.grid_res;
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

             ensure_grid();

             const TPMGrid& grid = grid_;
             const int nc = grid.n_combi;
             std::vector<double>& TPM_pred = TPM_pred_;

             const double   inv_nm1 = 1.0 / (nc - 1.0);

             if (row_start == 0) {
               TPM_post_pending_.assign(nc, 0.0);
               first_active_  = true;
               prev_active_r_ = -1;
               init_streams(n);
               post_history_.assign(n, std::vector<double>{});
               prev_active_r_history_.assign(n, -1);
             }

             for(int r = row_start; r < end; ++r) {
               if (at_mask[r]) {
                 post_history_[r] = TPM_post_pending_;
                 prev_active_r_history_[r] = prev_active_r_;

                 const double cp      = cp_col[r];
                 const double x       = cov_ptr[r];
                 const bool   reset   = first_active_ || (belief_reset_ && belief_reset_[r]);
                 const bool   prev_na = (prev_active_r_ < 0) || !nm[prev_active_r_];
                 const int    prev    = prev_na ? -1 : static_cast<int>(cov_ptr[prev_active_r_]);
                 const bool   curr_na = !nm[r];
                 const int    curr    = curr_na ? -1 : static_cast<int>(x);


                 if(reset) {
                   ensure_prior(a0_col[r], b0_col[r]);
                   TPM_pred = TPM_prior_;
                 } else {
                   // this is all correct, but the loops over the grid can be vectorised:
                   // const double sum_post = std::accumulate(TPM_post_pending_.begin(), TPM_post_pending_.end(), 0.0);
                   // for(int k = 0; k < nc; ++k) TPM_pred[k] = (1.0 - cp) * TPM_post_pending_[k] + cp * (sum_post - TPM_post_pending_[k]) * inv_nm1;
                   // normalise_inplace(TPM_pred);  // prior normalisation is done in ensure_prior()

                   // vectorised std::accumulate:
                   double sum_post = 0.0;
#pragma omp simd reduction(+:sum_post)
                   for (int k = 0; k < nc; ++k) {
                     sum_post += TPM_post_pending_[k];
                   }

                   // vectorised grid update:
                   const double keep = 1.0 - cp;
#pragma omp simd
                   for(int k = 0; k < nc; ++k) {
                     const double post = TPM_post_pending_[k];
                     TPM_pred[k] = keep * post + cp * (sum_post - post) * inv_nm1;
                   }
                   // and here we could also still pass a sum to get rid of one of the loops in normalise_inplace
                   normalise_inplace(TPM_pred);
                 }

                 double mean, mode, log_precision;
                 if (prev_na || reset) {
                   discrete_summary(grid.mean_p, TPM_pred, mean, mode, log_precision);  // calculated in-place; shares arithmic
                   set_active(r, mean, mode, log_precision);
                 } else {
                   const std::vector<double>& g = (prev == 1) ? grid.p_XX : grid.p_XY;
                   discrete_summary(g, TPM_pred, mean, mode, log_precision);  // calculated in-place; shares arithmic
                   set_active(r, mean, mode, log_precision);
                 }
                 if (!curr_na) set_surprise(r, x);

                 // Prediction summaries are finished.
                 // Transfer prediction into the pending posterior without copying.
                 // The old posterior becomes prediction scratch space for the next row.
                 TPM_post_pending_.swap(TPM_pred);

                 // compute pending posterior
                 if (!curr_na && !prev_na && !reset) {
                   const std::vector<double>& like = (prev == 0) ? (curr == 0 ? grid.like_YY : grid.like_XY) : (curr == 0 ? grid.like_YX : grid.like_XX);

                   double sum_update = 0.0;
                   // One pass: likelihood multiplication + normalisation sum.
#pragma omp simd reduction(+:sum_update)
                   for (int k = 0; k < nc; ++k) {
                     const double updated = TPM_post_pending_[k] * like[k];
                     TPM_post_pending_[k] = updated;
                     sum_update += updated;
                   }

                   // Second pass: normalise.
                   if (sum_update > 0.0) {
                     const double inv = 1.0 / sum_update;
#pragma omp simd
                     for (int k = 0; k < nc; ++k) {
                       TPM_post_pending_[k] *= inv;
                     }
                   }
                 }

                 first_active_  = false;
                 prev_active_r_ = r;
               } else {
                 set_inactive(r);
               }
             }
             mark_run_complete(end);
           }
};


// ---- Type mapping + factory ----

KernelType to_kernel_type(const Rcpp::String& k);

std::unique_ptr<BaseKernel> make_kernel(KernelType kt,
                                        SEXP custom_fun = R_NilValue);


#endif // KERNELS_H

