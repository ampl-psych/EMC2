#include <unordered_map>
#include <memory>
#include <array>
#include <vector>     //
#include <Rcpp.h>    //
#include "EMC2/userfun.hpp"

// View
struct KernelParsView {
  int n_rows;
  std::vector<const double*> cols;  // cols[k][row] = value for param k at trial row
};

// ---- Types ----

enum class KernelType {
  SimpleDelta,
  Delta2Kernel,
  // Delta2Kernel2,
  Delta2LR,
  PierceHall,
  LinIncr,
  LinDecr,
  ExpIncr,
  ExpDecr,
  PowIncr,
  PowDecr,
  Poly2,
  Poly3,
  Poly4,
  Custom
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
  case KernelType::PierceHall:
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

public:
  virtual ~BaseKernel() {}

  virtual void run(const KernelParsView& kernel_pars,
                   const Rcpp::NumericMatrix& covariate,
                   const std::vector<int>& comp_idx) = 0;

  virtual void reset() {
    out_.clear();
    has_run_ = false;
    expand_idx_.clear();
    has_expand_idx_ = false;
  }

  bool has_run() const { return has_run_; }

  const std::vector<double>& get_output() const { return out_; }

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

  // single-stream getter, code=1 for main trajectory by default. code=2 for pes in delta, code=3 for xx in new kernels
  virtual Rcpp::NumericVector get_output_stream(int code) const {
    using namespace Rcpp;
    if (code != 1) {
      stop("BaseKernel::get_output_stream: unsupported code %d (only 1)", code);
    }
    // out_ is already full-length at this point
    return wrap(out_);  // copies to NumericVector
  }

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

public:
  virtual ~DeltaKernel() {}

  const std::vector<double>& get_pes() const {
    return pes_;
  }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 2);
  }

  Rcpp::NumericVector get_output_stream(int code) const override {
    using namespace Rcpp;

    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      // main trajectory: already full-length
      return wrap(out_);
    }

    if (code == 2) {
      NumericVector res(n_full);

      if (!has_expand_idx_) {
        // no 'at': one-to-one
        if ((int)pes_.size() != n_full)
          stop("DeltaKernel: pes_ length mismatch");
        for (int i = 0; i < n_full; ++i) res[i] = pes_[i];
      } else {
        // with 'at': expand from compressed index
        const auto& idx = expand_idx_;
        if ((int)idx.size() != n_full)
          stop("DeltaKernel: expand_idx length mismatch");
        for (int i = 0; i < n_full; ++i) {
          int k = idx[i] - 1;  // compressed index
          res[i] = pes_[k];
        }
      }
      return res;
    }

    stop("DeltaKernel::get_output_stream: unsupported code %d (1=Q,2=PE)", code);
  }

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

             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate(r,0);
               if (!ISNAN(x)) {
                 out_[j] = x;  // compressed index
               }
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

             // double last = NA_REAL;
             for (int j = 0; j < n_comp; ++j) {
               int r = comp_idx[j];
               double x = covariate(r,0);
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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
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
           const Rcpp::NumericMatrix& covariate,
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

               double x = covariate(r,0);
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
               double x = covariate(r,0);
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

//
struct PierceHall : DeltaKernel {
  // Based on https://pmc.ncbi.nlm.nih.gov/articles/PMC4563025/#FD2
  // alpha_t+1 = alpha_t*abs(PE_t)*eta + (1-eta)*alpha_t
  // Three parameters: q0, alpha0, and eta

  // Only works if PEs are in range (-1, 1)!
  std::vector<double> alphas_;
  double alpha_;

  PierceHall() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("PierceHall expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, NA_REAL);
             pes_.assign(n_comp, NA_REAL);
             alphas_.assign(n_comp, NA_REAL);

             const double* q0_col       = kernel_pars.cols[0];
             const double* alpha0_col   = kernel_pars.cols[1];
             const double* eta_col      = kernel_pars.cols[2];

             int row0 = comp_idx[0];
             out_[0] = q_ = q0_col[row0];
             alphas_[0] = alpha_ = alpha0_col[row0];

             double pe = NA_REAL;

             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];
               double x = covariate(r,0);
               if (!ISNAN(x)) {
                 pe = x - q_;
                 q_ += alpha_ * pe;

                 // update alpha for next trial
                 double eta = eta_col[r];
                 alpha_ = std::abs(pe) * eta + (1-eta) * alpha_;
               } else {
                 pe = NA_REAL;
               }

               pes_[j] = pe;           // compressed index
               alphas_[j+1] = alpha_;
               out_[j+1] = q_;
             }

             mark_run_complete();
           }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 3);
  }

  Rcpp::NumericVector get_output_stream(int code) const override {
    using namespace Rcpp;

    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      // main trajectory: already full-length
      return wrap(out_);
    }

    // For other codes, choose the compressed source vector
    const std::vector<double>* src = nullptr;
    if (code == 2) {            // PE
      src = &pes_;
    } else if (code == 3) {     // learning rates
      src = &alphas_;
    } else {
      stop("PierceHall::get_output_stream: unsupported code %d "
             "(1=Q,2=PE,3=alpha)", code);
    }

    NumericVector res(n_full);

    if (!has_expand_idx_) {
      // no 'at': compressed and full coincide
      if ((int)src->size() != n_full) {
        stop("PierceHall::get_output_stream: source length (%d) != n_full (%d)",
             (int)src->size(), n_full);
      }
      for (int i = 0; i < n_full; ++i) {
        res[i] = (*src)[i];
      }
    } else {
      // with 'at': expand from compressed to full using expand_idx_
      const auto& idx = expand_idx_;
      if ((int)idx.size() != n_full) {
        stop("PierceHall::get_output_stream: expand_idx length (%d) != n_full (%d)",
             (int)idx.size(), n_full);
      }
      const int n_comp = static_cast<int>(src->size());
      for (int i = 0; i < n_full; ++i) {
        int k = idx[i] - 1;  // 1-based -> 0-based compressed index
        if (k < 0 || k >= n_comp) {
          stop("PierceHall::get_output_stream: index %d out of range [0,%d)", k, n_comp);
        }
        res[i] = (*src)[k];
      }

      return res;
    }

    stop("PierceHall::get_output_stream: unsupported code %d (1=Q,2=PEs,3=alpha)", code);
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qvalue";
    if (code == 2) return "PEs";
    if (code == 3) return "alpha";
    throw std::runtime_error("PierceHall::output_stream_name: unsupported code");
  }
};

// 2D PE kernel: separate from DeltaKernel
struct Delta2Kernel : SequentialKernel {
  double qFast_ = NA_REAL;
  double qSlow_ = NA_REAL;
  double q_     = NA_REAL;

  // [compressed trial][0 = fast PE, 1 = slow PE]
  std::vector<double> pes_fast_;
  std::vector<double> pes_slow_;
  std::vector<double> q_fast_;
  std::vector<double> q_slow_;

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
               double x = covariate(r,0);
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

  Rcpp::NumericVector get_output_stream(int code) const override {
    using namespace Rcpp;

    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      // main trajectory: already full-length
      return wrap(out_);
    }

    // For other codes, choose the compressed source vector
    const std::vector<double>* src = nullptr;
    if (code == 2) {            // Qfast
      src = &q_fast_;
    } else if (code == 3) {     // Qslow
      src = &q_slow_;
    } else if (code == 4) {     // PEfast
      src = &pes_fast_;
    } else if (code == 5) {     // PEslow
      src = &pes_slow_;
    } else {
      stop("Delta2Kernel::get_output_stream: unsupported code %d "
             "(1=Q,2=Qfast,3=Qslow,4=PEfast,5=PEslow)", code);
    }

    NumericVector res(n_full);

    if (!has_expand_idx_) {
      // no 'at': compressed and full coincide
      if ((int)src->size() != n_full) {
        stop("Delta2Kernel::get_output_stream: source length (%d) != n_full (%d)",
             (int)src->size(), n_full);
      }
      for (int i = 0; i < n_full; ++i) {
        res[i] = (*src)[i];
      }
    } else {
      // with 'at': expand from compressed to full using expand_idx_
      const auto& idx = expand_idx_;
      if ((int)idx.size() != n_full) {
        stop("Delta2Kernel::get_output_stream: expand_idx length (%d) != n_full (%d)",
             (int)idx.size(), n_full);
      }
      const int n_comp = static_cast<int>(src->size());
      for (int i = 0; i < n_full; ++i) {
        int k = idx[i] - 1;  // 1-based -> 0-based compressed index
        if (k < 0 || k >= n_comp) {
          stop("Delta2Kernel::get_output_stream: index %d out of range [0,%d)", k, n_comp);
        }
        res[i] = (*src)[k];
      }

      return res;
    }

    stop("Delta2Kernel::get_output_stream: unsupported code %d (1=Q,2=Qfast,3=Qslow,4=PEfast,5=PEslow)", code);
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


