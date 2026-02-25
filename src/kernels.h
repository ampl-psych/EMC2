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
  PearceHall,
  DeltaRisk,
  LinIncr,
  LinDecr,
  ExpIncr,
  ExpDecr,
  PowIncr,
  PowDecr,
  Poly2,
  Poly3,
  Poly4,
  VKFBinary,
  VKFContinous,
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
  case KernelType::PearceHall:
  case KernelType::DeltaRisk:
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
  case KernelType::VKFBinary:
  case KernelType::VKFContinous:
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

  // Step-wise updating methods
  // By default, kernels do not support external step-wise orchestration
  virtual bool supports_stepwise() const { return false; }

  // Called once per sequence before stepping.
  // Default: fallback to batch mode so non-stepwise kernels still work.
  virtual void init_stepwise(const KernelParsView& kernel_pars,
                             const Rcpp::NumericMatrix& covariate,
                             const std::vector<int>& comp_idx) {
    run(kernel_pars, covariate, comp_idx);
  }

  // Called once per compressed trial (j_comp) after init_stepwise().
  // Default: no-op.
  virtual void step(const KernelParsView& kernel_pars,
                    const Rcpp::NumericMatrix& covariate,
                    int r_full,   // full row index = comp_idx[j_comp]
                    int j_comp) {
    (void) kernel_pars;
    (void) covariate;
    (void) r_full;
    (void) j_comp;
  }

  //
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

             const double* cov_col0 = &covariate(0, 0);
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

               const double x = cov_col0[r];   // direct pointer access - faster
               if (!ISNAN(x)) {
                 double alpha = alpha_col[r];
                 pe = x - q_;
                 q_ += alpha * pe;
               } else {
                 pe = NA_REAL;
               }
               pes_[j] = pe;                   // compressed index
               out_[j+1] = q_;
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
// struct PearceHall : DeltaKernel {
//   // Based on https://pmc.ncbi.nlm.nih.gov/articles/PMC4563025/#FD2
//   // alpha_t+1 = alpha_t*abs(PE_t)*eta + (1-eta)*alpha_t
//   // Three parameters: q0, alpha0, and eta
//
//   // Only works if PEs are in range (-1, 1)!
//   std::vector<double> alphas_;
//   double alpha_;
//
//   PearceHall() {}
//
//   void run(const KernelParsView& kernel_pars,
//            const Rcpp::NumericMatrix& covariate,
//            const std::vector<int>& comp_idx) override {
//              if (kernel_pars.cols.size() != 3) {
//                Rcpp::stop("PearceHall expects 3 parameter columns, got %d",
//                           (int)kernel_pars.cols.size());
//              }
//
//              int n_comp = comp_idx.size();
//              out_.assign(n_comp, NA_REAL);
//              pes_.assign(n_comp, NA_REAL);
//              alphas_.assign(n_comp, NA_REAL);
//
//              const double* q0_col       = kernel_pars.cols[0];
//              const double* alpha0_col   = kernel_pars.cols[1];
//              const double* eta_col      = kernel_pars.cols[2];
//
//              int row0 = comp_idx[0];
//              out_[0] = q_ = q0_col[row0];
//              alphas_[0] = alpha_ = alpha0_col[row0];
//
//              double pe = NA_REAL;
//
//              for (int j = 0; j < n_comp - 1; ++j) {
//                int r = comp_idx[j];
//                double x = covariate(r,0);
//                if (!ISNAN(x)) {
//                  pe = x - q_;
//                  q_ += alpha_ * pe;
//
//                  // update alpha for next trial
//                  double eta = eta_col[r];
//                  alpha_ = std::abs(pe) * eta + (1-eta) * alpha_;
//                } else {
//                  pe = NA_REAL;
//                }
//
//                pes_[j] = pe;           // compressed index
//                alphas_[j+1] = alpha_;
//                out_[j+1] = q_;
//              }
//
//              mark_run_complete();
//            }
//
//   bool has_output_stream(int code) const override {
//     return (code >= 1 && code <= 3);
//   }
//
//   Rcpp::NumericVector get_output_stream(int code) const override {
//     using namespace Rcpp;
//
//     const int n_full = static_cast<int>(out_.size());
//
//     if (code == 1) {
//       // main trajectory: already full-length
//       return wrap(out_);
//     }
//
//     // For other codes, choose the compressed source vector
//     const std::vector<double>* src = nullptr;
//     if (code == 2) {            // PE
//       src = &pes_;
//     } else if (code == 3) {     // learning rates
//       src = &alphas_;
//     } else {
//       stop("PearceHall::get_output_stream: unsupported code %d "
//              "(1=Q,2=PE,3=alpha)", code);
//     }
//
//     NumericVector res(n_full);
//
//     if (!has_expand_idx_) {
//       // no 'at': compressed and full coincide
//       if ((int)src->size() != n_full) {
//         stop("PearceHall::get_output_stream: source length (%d) != n_full (%d)",
//              (int)src->size(), n_full);
//       }
//       for (int i = 0; i < n_full; ++i) {
//         res[i] = (*src)[i];
//       }
//     } else {
//       // with 'at': expand from compressed to full using expand_idx_
//       const auto& idx = expand_idx_;
//       if ((int)idx.size() != n_full) {
//         stop("PearceHall::get_output_stream: expand_idx length (%d) != n_full (%d)",
//              (int)idx.size(), n_full);
//       }
//       const int n_comp = static_cast<int>(src->size());
//       for (int i = 0; i < n_full; ++i) {
//         int k = idx[i] - 1;  // 1-based -> 0-based compressed index
//         if (k < 0 || k >= n_comp) {
//           stop("PearceHall::get_output_stream: index %d out of range [0,%d)", k, n_comp);
//         }
//         res[i] = (*src)[k];
//       }
//
//       return res;
//     }
//
//     stop("PearceHall::get_output_stream: unsupported code %d (1=Q,2=PEs,3=alpha)", code);
//   }
//
//   std::string output_stream_name(int code) const override {
//     if (code == 1) return "Qvalue";
//     if (code == 2) return "PEs";
//     if (code == 3) return "alpha";
//     throw std::runtime_error("PearceHall::output_stream_name: unsupported code");
//   }
// };


struct PearceHall : DeltaKernel {
  // Based on https://pmc.ncbi.nlm.nih.gov/articles/PMC4563025/#FD2
  // alpha_{t+1} = alpha_t * abs(PE_t) * eta + (1 - eta) * alpha_t
  // Parameters: q0, alpha0, eta

  std::vector<double> alphas_;  // per-trial learning rates
  double alpha_;                // current learning rate

  // step-wise / shared-latent support
  std::vector<int> comp_idx_;   // compressed index -> full trial
  int n_comp_ = 0;

  bool use_shared_alpha_ = false;
  std::vector<double>* shared_alpha_ = nullptr; // owned by TrendOpRuntime

  PearceHall() {}

  // Called from TrendRuntime::bind_all_ops_to_paramtable when shared_latent=TRUE
  void enable_shared_alpha(std::vector<double>* shared_alpha_vec) {
    use_shared_alpha_ = (shared_alpha_vec != nullptr);
    shared_alpha_ = shared_alpha_vec;
  }

  bool supports_stepwise() const override {
    // We only need external stepwise orchestration when we share alpha.
    return use_shared_alpha_;
  }

  void reset() override {
    DeltaKernel::reset();
    alphas_.clear();
    alpha_ = NA_REAL;
    comp_idx_.clear();
    n_comp_ = 0;
    // NOTE: shared_alpha_ is owned by TrendOpRuntime; do not clear here.
  }

  // -------- Original batch run (independent mode) --------
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             if (kernel_pars.cols.size() != 3) {
               Rcpp::stop("PearceHall expects 3 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, NA_REAL);
             pes_.assign(n_comp, NA_REAL);
             alphas_.assign(n_comp, NA_REAL);

             const double* q0_col     = kernel_pars.cols[0];
             const double* alpha0_col = kernel_pars.cols[1];
             const double* eta_col    = kernel_pars.cols[2];

             int row0 = comp_idx[0];
             out_[0]    = q_     = q0_col[row0];
             alphas_[0] = alpha_ = alpha0_col[row0];

             double pe = NA_REAL;

             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];
               double x = covariate(r, 0);
               if (!ISNAN(x)) {
                 pe = x - q_;
                 q_ += alpha_ * pe;

                 double eta = eta_col[r];
                 alpha_ = std::abs(pe) * eta + (1 - eta) * alpha_;
               } else {
                 pe = NA_REAL;
               }

               pes_[j]      = pe;      // compressed index
               alphas_[j+1] = alpha_;
               out_[j+1]    = q_;
             }

             mark_run_complete();
           }

  // -------- Step-wise API (used when shared_latent=TRUE) --------
  void init_stepwise(const KernelParsView& kernel_pars,
                     const Rcpp::NumericMatrix& covariate,
                     const std::vector<int>& comp_idx) override {

       if (kernel_pars.cols.size() != 3) {
         Rcpp::stop("PearceHall expects 3 parameter columns, got %d",
                    (int)kernel_pars.cols.size());
       }

       comp_idx_ = comp_idx;
       n_comp_   = static_cast<int>(comp_idx_.size());

       if (n_comp_ <= 0) {
         // nothing to run, clear & set complete
         out_.clear();
         pes_.clear();
         alphas_.clear();
         if (shared_alpha_) shared_alpha_->clear();
         mark_run_complete();
         return;
       }

       // Initialize
       out_.assign(n_comp_, NA_REAL);
       pes_.assign(n_comp_, NA_REAL);
       alphas_.assign(n_comp_, NA_REAL); // note that this one is shared

       const double* q0_col     = kernel_pars.cols[0];
       const double* alpha0_col = kernel_pars.cols[1];

       int row0   = comp_idx_[0];
       q_         = q0_col[row0];
       alpha_     = alpha0_col[row0];

       out_[0]    = q_;
       alphas_[0] = alpha_;

       // Initialize shared alpha trace (compressed index) if provided
       if (shared_alpha_) {
         if ((int)shared_alpha_->size() != n_comp_) {
           shared_alpha_->assign(n_comp_, NA_REAL);
         }
         (*shared_alpha_)[0] = alpha_;
       }

       // Note: we don't iterate here; actual updates happen in step()
       has_run_ = false; // will be set true after last step
     }

  void step(const KernelParsView& kernel_pars,
            const Rcpp::NumericMatrix& covariate,
            int r_full,
            int j_comp) override {

              using namespace Rcpp;

              if (!use_shared_alpha_) {
                // In non-shared mode we don't expect external stepwise calls.
                return;
              }

              if (n_comp_ <= 0) return;

              // j_comp is current compressed index; we write to j_comp+1 if possible.
              if (j_comp >= n_comp_ - 1) {
                // nothing to write to (no j_comp+1); mark done on last call
                if (!has_run_) mark_run_complete();
                return;
              }

              const double* eta_col = kernel_pars.cols[2];
              const int next = j_comp + 1;

              // 1) Get alpha_t from shared_alpha_ if available, else keep local alpha_
              double alpha_t = alpha_;
              if (shared_alpha_ && j_comp < (int)shared_alpha_->size()) {
                double a = (*shared_alpha_)[j_comp];
                if (!ISNAN(a)) alpha_t = a;
              }
              alpha_ = alpha_t;

              // 2) Read outcome
              double x  = covariate(r_full, 0);
              double pe = NA_REAL;

              if (!ISNAN(x)) {
                // ---- ACTIVE (chosen) covariate: do learning and write shared alpha ----
                pe = x - q_;                       // PE_t
                q_ += alpha_ * pe;                 // Q update

                double eta = eta_col[r_full];
                alpha_ = std::abs(pe) * eta + (1 - eta) * alpha_;   // alpha_{t+1}

                // store trajectories at t+1
                pes_[j_comp]  = pe;
                alphas_[next] = alpha_;
                out_[next]    = q_;

                if (shared_alpha_ && next < (int)shared_alpha_->size()) {
                  (*shared_alpha_)[next] = alpha_;  // write shared alpha_{t+1}
                }

              } else {
                // ---- PASSIVE (unchosen) covariate: NO learning, NO shared-alpha write ----
                pes_[j_comp] = NA_REAL;

                // Try to mirror whatever shared alpha_{t+1} is (if already written),
                // otherwise carry forward alpha_t.
                double alpha_next = alpha_;
                if (shared_alpha_ && next < (int)shared_alpha_->size()) {
                  double a_next = (*shared_alpha_)[next];
                  if (!ISNAN(a_next)) {
                    alpha_next = a_next;   // another kernel has already updated alpha_{t+1}
                  }
                }
                alphas_[next] = alpha_next;
                out_[next]    = q_;        // Q stays constant for this covariate
              }

              if (next == n_comp_ - 1) {
                mark_run_complete();
              }
            }

  // -------- Output streams --------
  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 3);
  }

  Rcpp::NumericVector get_output_stream(int code) const override {
    using namespace Rcpp;

    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      return wrap(out_);
    }

    const std::vector<double>* src = nullptr;
    if (code == 2) {            // PE
      src = &pes_;
    } else if (code == 3) {     // learning rates
      src = &alphas_;
    } else {
      stop("PearceHall::get_output_stream: unsupported code %d "
             "(1=Q,2=PE,3=alpha)", code);
    }

    NumericVector res(n_full);

    if (!has_expand_idx_) {
      if ((int)src->size() != n_full) {
        stop("PearceHall::get_output_stream: source length (%d) != n_full (%d)",
             (int)src->size(), n_full);
      }
      for (int i = 0; i < n_full; ++i) {
        res[i] = (*src)[i];
      }
      return res;
    } else {
      const auto& idx = expand_idx_;
      if ((int)idx.size() != n_full) {
        stop("PearceHall::get_output_stream: expand_idx length (%d) != n_full (%d)",
             (int)idx.size(), n_full);
      }
      const int n_comp = static_cast<int>(src->size());
      for (int i = 0; i < n_full; ++i) {
        int k = idx[i] - 1;  // 1-based -> 0-based compressed index
        if (k < 0 || k >= n_comp) {
          stop("PearceHall::get_output_stream: index %d out of range [0,%d)", k, n_comp);
        }
        res[i] = (*src)[k];
      }
      return res;
    }
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qvalue";
    if (code == 2) return "PEs";
    if (code == 3) return "alpha";
    throw std::runtime_error("PearceHall::output_stream_name: unsupported code");
  }
};

struct DeltaRisk : DeltaKernel {
  // Based on d'Acremont et al, 2009: doi.org/10.1016/j.neuroimage.2009.04.096
  // PE_t = cov - Q_{t}
  // Q_t+1 = Q_{t} + alpha * PE/sqrt(h_t)

  // h_t is an estimate of "risk" (i.e., variance); learnt via delta rule learning:
  // Risk prediction error xi_t:
  // xi_t = (PE_t)^2 - h_t
  // h_t+t = h_{t} + alpha_risk * xi_t

  // Four parameters: q0, alpha, xi_0, alpha_risk

  std::vector<double> xis_;  // risk PEs
  double xi; // risk PE
  std::vector<double> hs_;  // risk estimates
  double h_; // risk estimate

  DeltaRisk() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {
             if (kernel_pars.cols.size() != 4) {
               Rcpp::stop("DeltaRisk expects 4 parameter columns, got %d",
                          (int)kernel_pars.cols.size());
             }

             int n_comp = comp_idx.size();
             out_.assign(n_comp, NA_REAL);
             pes_.assign(n_comp, NA_REAL);
             hs_.assign(n_comp, NA_REAL);   // risk estimates
             xis_.assign(n_comp, NA_REAL);  // risk prediction errors

             const double* q0_col       = kernel_pars.cols[0];
             const double* alpha_col    = kernel_pars.cols[1];
             const double* h0_col       = kernel_pars.cols[2];
             const double* alpha_risk_col = kernel_pars.cols[3];

             int row0 = comp_idx[0];
             out_[0] = q_ = q0_col[row0];
             hs_[0] = h_ = h0_col[row0];

             double pe = NA_REAL;

             for (int j = 0; j < n_comp - 1; ++j) {
               int r = comp_idx[j];
               double x = covariate(r,0);
               if (!ISNAN(x)) {
                 double alpha_ = alpha_col[r];
                 double alpha_risk = alpha_risk_col[r];

                 pe = x - q_;
                 q_ += alpha_ * pe / std::sqrt(h_);

                 // Risk prediction error
                 xi = (pe * pe) - h_;
                 h_ += alpha_risk * xi;

               } else {
                 pe = NA_REAL;
                 xi = NA_REAL;
               }

               pes_[j] = pe;           // compressed index
               xis_[j] = xi;           // compressed index
               hs_[j+1] = h_;            // compressed index
               out_[j+1] = q_;
             }

             mark_run_complete();
           }

  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 4);
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
    } else if (code == 3) {     // risk estimates
      src = &hs_;
    } else if (code == 4) {     // risk PEs
      src = &xis_;
    } else {
      stop("DeltaRisk::get_output_stream: unsupported code %d "
             "(1=Q,2=PE,3=alpha)", code);
    }

    NumericVector res(n_full);

    if (!has_expand_idx_) {
      // no 'at': compressed and full coincide
      if ((int)src->size() != n_full) {
        stop("DeltaRisk::get_output_stream: source length (%d) != n_full (%d)",
             (int)src->size(), n_full);
      }
      for (int i = 0; i < n_full; ++i) {
        res[i] = (*src)[i];
      }
    } else {
      // with 'at': expand from compressed to full using expand_idx_
      const auto& idx = expand_idx_;
      if ((int)idx.size() != n_full) {
        stop("DeltaRisk::get_output_stream: expand_idx length (%d) != n_full (%d)",
             (int)idx.size(), n_full);
      }
      const int n_comp = static_cast<int>(src->size());
      for (int i = 0; i < n_full; ++i) {
        int k = idx[i] - 1;  // 1-based -> 0-based compressed index
        if (k < 0 || k >= n_comp) {
          stop("DeltaRisk::get_output_stream: index %d out of range [0,%d)", k, n_comp);
        }
        res[i] = (*src)[k];
      }

      return res;
    }

    stop("DeltaRisk::get_output_stream: unsupported code %d (1=Q,2=PEs,3=risk estimate,4=risk PE)", code);
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "Qvalue";
    if (code == 2) return "PEs";
    if (code == 3) return "riskvalue";
    if (code == 4) return "riskPE";
    throw std::runtime_error("DeltaRisk::output_stream_name: unsupported code");
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

struct VKFKernel : SequentialKernel {
  // Latent state at latest compressed step
  double m_      = NA_REAL;  // prediction (latent state, linear scale)
  double v_      = NA_REAL;  // volatility
  double w_      = NA_REAL;  // uncertainty (posterior variance)
  double sigma2_ = NA_REAL;  // observation noise variance

  // Per-trial (compressed) trajectories
  std::vector<double> ms_;     // m_t
  std::vector<double> vs_;     // v_t
  std::vector<double> ws_;     // w_t
  std::vector<double> pes_m_;  // delta_m
  std::vector<double> pes_v_;  // delta_v
  std::vector<double> ks_;     // k_t (Kalman gain / learning rate)

  VKFKernel() {}

  void reset() override {
    SequentialKernel::reset();
    m_ = v_ = w_ = sigma2_ = NA_REAL;
    ms_.clear();
    vs_.clear();
    ws_.clear();
    pes_m_.clear();
    pes_v_.clear();
    ks_.clear();
  }

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             using namespace Rcpp;

             // Expect: m0, v0, w0, lambda, sigma2
             if (kernel_pars.cols.size() != 5) {
               stop("VKFKernel expects 5 parameter columns (m0, v0, w0, lambda, sigma2), got %d",
                    (int)kernel_pars.cols.size());
             }

             const int n_comp = static_cast<int>(comp_idx.size());
             if (n_comp <= 0) {
               out_.clear();
               ms_.clear();
               vs_.clear();
               ws_.clear();
               pes_m_.clear();
               pes_v_.clear();
               ks_.clear();
               mark_run_complete();
               return;
             }

             // Allocate compressed trajectories
             out_.assign(n_comp, NA_REAL);   // main: m_t (prediction)
             ms_.assign(n_comp, NA_REAL);
             vs_.assign(n_comp, NA_REAL);
             ws_.assign(n_comp, NA_REAL);
             pes_m_.assign(n_comp, NA_REAL);
             pes_v_.assign(n_comp, NA_REAL);
             ks_.assign(n_comp, NA_REAL);

             // Parameter columns
             const double* m0_col     = kernel_pars.cols[0];
             const double* v0_col     = kernel_pars.cols[1];
             const double* w0_col     = kernel_pars.cols[2];
             const double* lambda_col = kernel_pars.cols[3];
             const double* sigma2_col = kernel_pars.cols[4];

             // Initial compressed element (at comp_idx[0])
             const int row0 = comp_idx[0];
             m_      = m0_col[row0];
             v_      = v0_col[row0];
             w_      = w0_col[row0];
             sigma2_ = sigma2_col[row0];  // separate param, no w0 = sigma2 assumption here

             // Store initial state: index 0
             out_[0] = m_;
             ms_[0]  = m_;
             vs_[0]  = v_;
             ws_[0]  = w_;

             // Iterate over compressed trials j = 0..n_comp-2
             for (int j = 0; j < n_comp - 1; ++j) {
               const int r = comp_idx[j];
               const double o = covariate(r, 0);      // continuous outcome

               double pe_m = NA_REAL;
               double pe_v = NA_REAL;
               double k    = NA_REAL;

               if (!ISNAN(o)) {
                 const double lambda = lambda_col[r];
                 const double sigma2 = sigma2_col[r]; // can be trial-varying if needed

                 const double mpre = m_;
                 const double wpre = w_;

                 // Prediction error: o - m
                 pe_m = o - m_;

                 // Kalman gain / learning rate, Eq 9
                 const double num = w_ + v_;              // w_{t-1} + v_{t-1}
                 const double den = num + sigma2;         // + sigma^2
                 k = (den > 0.0) ? (num / den) : 0.0;

                 // Mean update, Eq 10
                 m_ = m_ + k * pe_m;

                 // Variance update, Eq 11
                 w_ = (1.0 - k) * (w_ + v_);

                 // Covariance term, Eq 12
                 const double wcov = (1.0 - k) * wpre;

                 // Volatility prediction error, Eq 13 inner term
                 pe_v = (m_ - mpre) * (m_ - mpre)
                   + w_ + wpre
                 - 2.0 * wcov
                 - v_;

                 // Volatility update, Eq 13
                 v_ = v_ + lambda * pe_v;
               }

               // Store trial-level values at compressed index j
               pes_m_[j] = pe_m;
               pes_v_[j] = pe_v;
               ks_[j]    = k;

               // Updated state at index j+1
               ms_[j + 1]  = m_;
               vs_[j + 1]  = v_;
               ws_[j + 1]  = w_;
               out_[j + 1] = m_;
             }

             mark_run_complete();
           }

  // Output streams:
  // 1 = m (prediction)
  // 2 = prediction error delta_m
  // 3 = learning rate k (Kalman gain)
  // 4 = volatility v
  // 5 = volatility prediction error delta_v
  // 6 = uncertainty w
  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 6);
  }

  Rcpp::NumericVector get_output_stream(int code) const override {
    using namespace Rcpp;

    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      // main trajectory: m_t
      return wrap(out_);
    }

    // Pick source vector
    const std::vector<double>* src = nullptr;
    if      (code == 2) src = &pes_m_;
    else if (code == 3) src = &ks_;
    else if (code == 4) src = &vs_;
    else if (code == 5) src = &pes_v_;
    else if (code == 6) src = &ws_;
    else {
      stop("VKFKernel::get_output_stream: unsupported code %d "
             "(1=m,2=PE_m,3=k,4=volatility,5=vol_PE,6=uncertainty)", code);
    }

    NumericVector res(n_full);

    if (!has_expand_idx_) {
      // no 'at': compressed == full
      if ((int)src->size() != n_full) {
        stop("VKFKernel::get_output_stream: source length (%d) != n_full (%d)",
             (int)src->size(), n_full);
      }
      for (int i = 0; i < n_full; ++i) {
        res[i] = (*src)[i];
      }
      return res;
    } else {
      // with 'at': expand using expand_idx_
      const auto& idx = expand_idx_;
      if ((int)idx.size() != n_full) {
        stop("VKFKernel::get_output_stream: expand_idx length (%d) != n_full (%d)",
             (int)idx.size(), n_full);
      }
      const int n_comp = static_cast<int>(src->size());
      for (int i = 0; i < n_full; ++i) {
        int k_idx = idx[i] - 1;  // 1-based -> 0-based
        if (k_idx < 0 || k_idx >= n_comp) {
          stop("VKFKernel::get_output_stream: index %d out of range [0,%d)",
               k_idx, n_comp);
        }
        res[i] = (*src)[k_idx];
      }
      return res;
    }
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "m";
    if (code == 2) return "PE_m";
    if (code == 3) return "k";
    if (code == 4) return "volatility";
    if (code == 5) return "vol_PE";
    if (code == 6) return "uncertainty";
    throw std::runtime_error("VKFKernel::output_stream_name: unsupported code");
  }
};


// struct VKFBinaryKernel : SequentialKernel {
//   double m_     = NA_REAL;    // latent mean
//   double v_     = NA_REAL;    // volatility
//   double w_     = NA_REAL;    // uncertainty
//   double omega_ = NA_REAL;    // noise parameter
//
//   std::vector<double> ms_;     // latent m_t
//   std::vector<double> vs_;
//   std::vector<double> ws_;
//   std::vector<double> pes_m_;
//   std::vector<double> pes_v_;
//   std::vector<double> lrs_;    // alpha_t
//   // out_ will hold p_t = sigmoid(m_t)
//
//   void reset() override {
//     SequentialKernel::reset();
//     m_ = v_ = w_ = omega_ = NA_REAL;
//     ms_.clear();
//     vs_.clear();
//     ws_.clear();
//     pes_m_.clear();
//     pes_v_.clear();
//     lrs_.clear();
//   }
//
//   void run(const KernelParsView& kernel_pars,
//            const Rcpp::NumericMatrix& covariate,
//            const std::vector<int>& comp_idx) override {
//
//              using namespace Rcpp;
//
//              // Expect: m0, v0, w0, lambda, omega
//              if (kernel_pars.cols.size() != 5) {
//                stop("VKFBinaryKernel expects 5 parameter columns (m0, v0, w0, lambda, omega), got %d",
//                     (int)kernel_pars.cols.size());
//              }
//
//              const int n_comp = static_cast<int>(comp_idx.size());
//              if (n_comp <= 0) {
//                out_.clear();
//                ms_.clear();
//                vs_.clear();
//                ws_.clear();
//                pes_m_.clear();
//                pes_v_.clear();
//                lrs_.clear();
//                mark_run_complete();
//                return;
//              }
//
//              out_.assign(n_comp, NA_REAL);   // main: p_t = s(m_t)
//              ms_.assign(n_comp, NA_REAL);
//              vs_.assign(n_comp, NA_REAL);
//              ws_.assign(n_comp, NA_REAL);
//              pes_m_.assign(n_comp, NA_REAL);
//              pes_v_.assign(n_comp, NA_REAL);
//              lrs_.assign(n_comp, NA_REAL);
//
//              // parameter columns
//              const double* m0_col     = kernel_pars.cols[0];
//              const double* v0_col     = kernel_pars.cols[1];
//              const double* w0_col     = kernel_pars.cols[2];
//              const double* lambda_col = kernel_pars.cols[3];
//              const double* omega_col  = kernel_pars.cols[4];
//
//              // initial compressed element
//              const int row0 = comp_idx[0];
//              m_     = m0_col[row0];
//              v_     = v0_col[row0];
//              w_     = w0_col[row0];
//              omega_ = omega_col[row0];    // now separate from w0
//
//              ms_[0]  = m_;
//              vs_[0]  = v_;
//              ws_[0]  = w_;
//              out_[0] = 1.0 / (1.0 + std::exp(-m_));  // p_0
//
//              for (int j = 0; j < n_comp - 1; ++j) {
//                int r = comp_idx[j];
//                double o = covariate(r, 0);
//
//                double pe_m = NA_REAL;
//                double pe_v = NA_REAL;
//                double alpha = NA_REAL;
//
//                if (!ISNAN(o)) {
//                  double lambda = lambda_col[r];
//                  double omega  = omega_col[r];  // allows trial-varying omega if needed
//
//                  double mpre = m_;
//                  double wpre = w_;
//
//                  // prediction error on Bernoulli prob scale: o - s(m)
//                  double p = 1.0 / (1.0 + std::exp(-m_));
//                  pe_m = o - p;
//
//                  // Kalman-like gain (Eq 14)
//                  double num = w_ + v_;
//                  double den = num + omega;
//                  double k   = (den > 0.0) ? (num / den) : 0.0;
//
//                  // learning rate alpha_t (Eq 15)
//                  double sum_wv = w_ + v_;
//                  alpha = (sum_wv > 0.0) ? std::sqrt(sum_wv) : 0.0;
//
//                  // mean update (Eq 16)
//                  m_ = m_ + alpha * pe_m;
//
//                  // variance update (Eq 17)
//                  w_ = (1.0 - k) * (w_ + v_);
//
//                  // covariance and volatility update (Eqs 18–19)
//                  double wcov = (1.0 - k) * wpre;
//                  pe_v = (m_ - mpre) * (m_ - mpre)
//                    + w_ + wpre
//                  - 2.0 * wcov
//                  - v_;
//                  v_ = v_ + lambda * pe_v;
//                }
//
//                pes_m_[j] = pe_m;
//                pes_v_[j] = pe_v;
//                lrs_[j]   = alpha;
//
//                ms_[j + 1]  = m_;
//                vs_[j + 1]  = v_;
//                ws_[j + 1]  = w_;
//                out_[j + 1] = 1.0 / (1.0 + std::exp(-m_));  // p_t
//              }
//
//              mark_run_complete();
//            }
//
//   // streams: 1 = p, 2 = PE_m, 3 = alpha, 4 = v, 5 = vol_PE, 6 = w, 7 = m
//   bool has_output_stream(int code) const override {
//     return (code >= 1 && code <= 7);
//   }
//
//   Rcpp::NumericVector get_output_stream(int code) const override {
//     using namespace Rcpp;
//     const int n_full = static_cast<int>(out_.size());
//
//     if (code == 1) return wrap(out_);  // p_t
//
//     const std::vector<double>* src = nullptr;
//     if      (code == 2) src = &pes_m_;
//     else if (code == 3) src = &lrs_;
//     else if (code == 4) src = &vs_;
//     else if (code == 5) src = &pes_v_;
//     else if (code == 6) src = &ws_;
//     else if (code == 7) src = &ms_;
//     else stop("VKFBinaryKernel::get_output_stream: unsupported code");
//
//     Rcpp::NumericVector res(n_full);
//
//     if (!has_expand_idx_) {
//       if ((int)src->size() != n_full)
//         stop("VKFBinaryKernel::get_output_stream: source length mismatch");
//       for (int i = 0; i < n_full; ++i) res[i] = (*src)[i];
//       return res;
//     } else {
//       const auto& idx = expand_idx_;
//       if ((int)idx.size() != n_full)
//         stop("VKFBinaryKernel::get_output_stream: expand_idx length mismatch");
//       const int n_comp = static_cast<int>(src->size());
//       for (int i = 0; i < n_full; ++i) {
//         int k = idx[i] - 1;
//         if (k < 0 || k >= n_comp)
//           stop("VKFBinaryKernel::get_output_stream: index out of range");
//         res[i] = (*src)[k];
//       }
//       return res;
//     }
//   }
//
//   std::string output_stream_name(int code) const override {
//     if (code == 1) return "p";             // s(m)
//     if (code == 2) return "PE_m";
//     if (code == 3) return "alpha";
//     if (code == 4) return "volatility";
//     if (code == 5) return "vol_PE";
//     if (code == 6) return "uncertainty";
//     if (code == 7) return "m_latent";
//     throw std::runtime_error("VKFBinaryKernel::output_stream_name: unsupported code");
//   }
// };

// Binary VKF kernel with parameters: m0, alpha0, r, lambda, omega
// - m0     : initial latent state (not on probability scale; p = sigmoid(m0))
// - alpha0 : initial learning rate (alpha_0)
// - r      : fraction of initial variance attributed to volatility (0<r<1)
// - lambda : volatility learning rate
// - omega  : noise parameter for inference (acts like observation noise in k_t)
//
// Internally:
//   S0   = alpha0^2
//   v0   = r * S0
//   w0   = (1 - r) * S0
//
// Learning rate each trial: alpha_t = sqrt(w_t-1 + v_t-1)
// Kalman-like gain:         k_t     = (w_t-1 + v_t-1) / (w_t-1 + v_t-1 + omega)

struct VKFBinaryKernel : SequentialKernel {
  // Latent state at latest compressed step
  double m_     = NA_REAL;  // latent mean
  double v_     = NA_REAL;  // volatility
  double w_     = NA_REAL;  // uncertainty (posterior variance)
  double omega_ = NA_REAL;  // noise parameter

  // Per-trial (compressed) trajectories
  std::vector<double> ms_;     // m_t (latent)
  std::vector<double> vs_;     // v_t
  std::vector<double> ws_;     // w_t
  std::vector<double> pes_m_;  // delta_m
  std::vector<double> pes_v_;  // delta_v
  std::vector<double> lrs_;    // alpha_t

  // out_ (inherited) will hold p_t = sigmoid(m_t)

  VKFBinaryKernel() {}

  void reset() override {
    SequentialKernel::reset();
    m_ = v_ = w_ = omega_ = NA_REAL;
    ms_.clear();
    vs_.clear();
    ws_.clear();
    pes_m_.clear();
    pes_v_.clear();
    lrs_.clear();
  }

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericMatrix& covariate,
           const std::vector<int>& comp_idx) override {

             using namespace Rcpp;

             // Expect parameter columns: m0, alpha0, r, lambda, omega
             if (kernel_pars.cols.size() != 5) {
               stop("BinaryVKFKernel expects 5 parameter columns "
                      "(m0, alpha0, r, lambda, omega), got %d",
                      (int)kernel_pars.cols.size());
             }

             const int n_comp = static_cast<int>(comp_idx.size());
             if (n_comp <= 0) {
               out_.clear();
               ms_.clear();
               vs_.clear();
               ws_.clear();
               pes_m_.clear();
               pes_v_.clear();
               lrs_.clear();
               mark_run_complete();
               return;
             }

             // Allocate compressed trajectories
             out_.assign(n_comp, NA_REAL);   // main: p_t = sigmoid(m_t)
             ms_.assign(n_comp, NA_REAL);
             vs_.assign(n_comp, NA_REAL);
             ws_.assign(n_comp, NA_REAL);
             pes_m_.assign(n_comp, NA_REAL);
             pes_v_.assign(n_comp, NA_REAL);
             lrs_.assign(n_comp, NA_REAL);

             // Parameter columns
             const double* m0_col     = kernel_pars.cols[0];
             const double* alpha0_col = kernel_pars.cols[1];
             const double* r_col      = kernel_pars.cols[2];
             const double* lambda_col = kernel_pars.cols[3];
             const double* omega_col  = kernel_pars.cols[4];

             // Initial compressed element (row0)
             const int row0 = comp_idx[0];

             const double m0     = m0_col[row0];
             const double alpha0 = alpha0_col[row0];
             const double r      = r_col[row0];
             const double omega0 = omega_col[row0];

             // Map (alpha0, r) -> (w0, v0)
             const double S0 = alpha0 * alpha0;      // S0 = alpha0^2
             const double v0 = r       * S0;
             const double w0 = (1.0-r) * S0;

             m_     = m0;
             v_     = v0;
             w_     = w0;
             omega_ = omega0;

             // Store initial state at compressed index 0
             ms_[0]  = m_;
             vs_[0]  = v_;
             ws_[0]  = w_;
             out_[0] = 1.0 / (1.0 + std::exp(-m_));   // p_0 = sigmoid(m_0)

             // Iterate over compressed trials j = 0..n_comp-2
             for (int j = 0; j < n_comp - 1; ++j) {
               const int r_idx = comp_idx[j];
               const double o  = covariate(r_idx, 0);  // binary outcome (0/1, may be NA)

               double pe_m   = NA_REAL;
               double pe_v   = NA_REAL;
               double alpha_t = NA_REAL;

               if (!ISNAN(o)) {
                 const double lambda_t = lambda_col[r_idx];
                 const double omega_t  = omega_col[r_idx];  // allow trial-varying omega if desired

                 const double mpre = m_;
                 const double wpre = w_;

                 // Bernoulli prediction probability and error
                 const double p = 1.0 / (1.0 + std::exp(-m_));
                 pe_m = o - p;

                 // Kalman-like gain, Eq 14: k_t = (w+v)/(w+v+omega)
                 const double sum_wv = w_ + v_;
                 const double den    = sum_wv + omega_t;
                 const double k      = (den > 0.0) ? (sum_wv / den) : 0.0;

                 // Learning rate alpha_t, Eq 15
                 alpha_t = (sum_wv > 0.0) ? std::sqrt(sum_wv) : 0.0;

                 // Mean update, Eq 16: m_t = m_{t-1} + alpha_t * (o_t - s(m_{t-1}))
                 m_ = m_ + alpha_t * pe_m;

                 // Variance update, Eq 17: w_t = (1 - k_t) * (w_{t-1} + v_{t-1})
                 w_ = (1.0 - k) * (w_ + v_);

                 // Covariance and volatility update, Eqs 18–19
                 const double wcov = (1.0 - k) * wpre;  // Eq 18

                 // Volatility prediction error, Eq 19 inner term
                 pe_v = (m_ - mpre) * (m_ - mpre)
                   + w_ + wpre
                 - 2.0 * wcov
                 - v_;

                 // Volatility update, Eq 19: v_t = v_{t-1} + lambda * delta_v
                 v_ = v_ + lambda_t * pe_v;
               }

               // Store signals at compressed index j
               pes_m_[j] = pe_m;
               pes_v_[j] = pe_v;
               lrs_[j]   = alpha_t;

               // Updated state at compressed index j+1
               ms_[j + 1]  = m_;
               vs_[j + 1]  = v_;
               ws_[j + 1]  = w_;
               out_[j + 1] = 1.0 / (1.0 + std::exp(-m_));  // p_t = sigmoid(m_t)
             }

             mark_run_complete();
           }

  // Output streams:
  //  1 = p (probability, sigmoid(m))
  //  2 = prediction error delta_m
  //  3 = alpha_t (learning rate)
  //  4 = volatility v
  //  5 = volatility prediction error delta_v
  //  6 = uncertainty w
  //  7 = m (latent state)
  bool has_output_stream(int code) const override {
    return (code >= 1 && code <= 7);
  }

  Rcpp::NumericVector get_output_stream(int code) const override {
    using namespace Rcpp;
    const int n_full = static_cast<int>(out_.size());

    if (code == 1) {
      // main: p_t = sigmoid(m_t)
      return wrap(out_);
    }

    const std::vector<double>* src = nullptr;
    if      (code == 2) src = &pes_m_;
    else if (code == 3) src = &lrs_;
    else if (code == 4) src = &vs_;
    else if (code == 5) src = &pes_v_;
    else if (code == 6) src = &ws_;
    else if (code == 7) src = &ms_;
    else {
      stop("BinaryVKFKernel::get_output_stream: unsupported code %d "
             "(1=p,2=PE_m,3=alpha,4=volatility,5=vol_PE,6=uncertainty,7=m)", code);
    }

    NumericVector res(n_full);

    if (!has_expand_idx_) {
      if ((int)src->size() != n_full) {
        stop("BinaryVKFKernel::get_output_stream: source length (%d) != n_full (%d)",
             (int)src->size(), n_full);
      }
      for (int i = 0; i < n_full; ++i) {
        res[i] = (*src)[i];
      }
      return res;
    } else {
      const auto& idx = expand_idx_;
      if ((int)idx.size() != n_full) {
        stop("BinaryVKFKernel::get_output_stream: expand_idx length (%d) != n_full (%d)",
             (int)idx.size(), n_full);
      }
      const int n_comp = static_cast<int>(src->size());
      for (int i = 0; i < n_full; ++i) {
        int k = idx[i] - 1;  // 1-based -> 0-based
        if (k < 0 || k >= n_comp) {
          stop("BinaryVKFKernel::get_output_stream: index %d out of range [0,%d)",
               k, n_comp);
        }
        res[i] = (*src)[k];
      }
      return res;
    }
  }

  std::string output_stream_name(int code) const override {
    if (code == 1) return "p";             // s(m)
    if (code == 2) return "PE_m";
    if (code == 3) return "alpha";
    if (code == 4) return "volatility";
    if (code == 5) return "vol_PE";
    if (code == 6) return "uncertainty";
    if (code == 7) return "m_latent";
    throw std::runtime_error("BinaryVKFKernel::output_stream_name: unsupported code");
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


