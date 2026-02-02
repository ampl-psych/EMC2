#include <unordered_map>
#include <memory>
#include <array>
#include <vector>     //
#include <Rcpp.h>    //

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
};

// ---- Base + hierarchy ----

struct BaseKernel {
protected:
  std::vector<double> out_;    // last run's trajectory
  bool has_run_ = false;

public:
  virtual ~BaseKernel() {}

  // kp: kernel_pars (compressed, n_comp x K)
  // covariate: compressed covariate vector (length n_comp)
  // comp_idx: 0-based indices into kp/covariate defining the block
  virtual void run(const KernelParsView& kernel_pars,
                   const Rcpp::NumericVector& covariate,
                   const Rcpp::IntegerVector& comp_idx) = 0;

  virtual void reset() {
    out_.clear();
    has_run_ = false;
  }

  bool has_run() const {
    return has_run_;
  }

  const std::vector<double>& get_output() const {
    return out_;
  }

protected:
  void mark_run_complete() { has_run_ = true; }
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
};


// ---- Individual kernels ----
// ---- Non-sequential kernels ----

struct LinIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);  // zeros, as in original code

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 out_[j] = cov_i;
               }
             }

             mark_run_complete();
           }
};

struct LinDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 out_[j] = -cov_i;
               }
             }

             mark_run_complete();
           }
};

struct ExpDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);

             const double* lambda_col = kernel_pars.cols[0];
             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 double lambda = lambda_col[row];
                 out_[j] = std::exp(-lambda * cov_i);
               }
             }

             mark_run_complete();
           }
};

struct ExpIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);
             const double* lambda_col = kernel_pars.cols[0];

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 double lambda = lambda_col[row];
                 out_[j] = 1.0 - std::exp(-lambda * cov_i);
               }
             }

             mark_run_complete();
           }
};

struct PowDecrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);

             const double* alpha_col = kernel_pars.cols[0];

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 double alpha = alpha_col[row];
                 out_[j] = std::pow(1.0 + cov_i, -alpha);
               }
             }

             mark_run_complete();
           }
};

struct PowIncrKernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);

             const double* alpha_col = kernel_pars.cols[0];

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 double alpha = alpha_col[row];
                 out_[j] = 1.0 - std::pow(1.0 + cov_i, -alpha);
               }
             }

             mark_run_complete();
           }
};

struct Poly2Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);

             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 double a1 = a1_col[row];
                 double a2 = a2_col[row];
                 out_[j] = a1 * cov_i + a2 * std::pow(cov_i, 2.0);
               }
             }

             mark_run_complete();
           }
};

struct Poly3Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);

             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];
             const double* a3_col = kernel_pars.cols[2];

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 double a1 = a1_col[row];
                 double a2 = a2_col[row];
                 double a3 = a3_col[row];
                 double cov2 = cov_i * cov_i;
                 out_[j] = a1 * cov_i + a2 * cov2 + a3 * cov2 * cov_i; // cov^3
               }
             }

             mark_run_complete();
           }
};

struct Poly4Kernel : BaseKernel {
  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {
             int n = comp_idx.size();
             out_.assign(n, 0.0);

             const double* a1_col = kernel_pars.cols[0];
             const double* a2_col = kernel_pars.cols[1];
             const double* a3_col = kernel_pars.cols[2];
             const double* a4_col = kernel_pars.cols[3];

             for (int j = 0; j < n; ++j) {
               int row = comp_idx[j];
               double cov_i = covariate[row];
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 double a1 = a1_col[row];
                 double a2 = a2_col[row];
                 double a3 = a3_col[row];
                 double a4 = a4_col[row];

                 double cov2 = cov_i * cov_i;
                 double cov3 = cov2 * cov_i;
                 double cov4 = cov2 * cov2;
                 out_[j] = a1 * cov_i + a2 * cov2 + a3 * cov3 + a4 * cov4;
               }
             }

             mark_run_complete();
           }
};


// Sequential kernels
struct SimpleDelta : DeltaKernel {
  SimpleDelta() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {

             int n_trials = comp_idx.size();

             out_.assign(n_trials, NA_REAL);
             pes_.assign(n_trials, NA_REAL);

             const double* q0_col = kernel_pars.cols[0];
             const double* alpha_col = kernel_pars.cols[1];


             int row0 = comp_idx[0];
             double q0 = q0_col[row0];
             q_ = q0;

             for (int j = 0; j < n_trials; ++j) {
               int row      = comp_idx[j];
               double cov_i = covariate[row];
               double alpha = alpha_col[row];

               double pe = NA_REAL;
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 pe = cov_i - q_;
                 q_ += alpha * pe;
               }

               pes_[j] = pe;
               out_[j] = q_;
             }

             mark_run_complete();
           }
};

struct Delta2LR : DeltaKernel {
  Delta2LR() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {

             int n_trials = comp_idx.size();

             out_.assign(n_trials, NA_REAL);
             pes_.assign(n_trials, NA_REAL);

             const double* q0_col = kernel_pars.cols[0];
             const double* alphaPos_col = kernel_pars.cols[1];
             const double* alphaNeg_col = kernel_pars.cols[2];

             int row0 = comp_idx[0];
             double q0 = q0_col[row0];
             q_ = q0;

             for (int j = 0; j < n_trials; ++j) {
               int row          = comp_idx[j];
               double cov_i     = covariate[row];
               double alphaPos  = alphaPos_col[row];
               double alphaNeg  = alphaNeg_col[row];

               double pe = NA_REAL;
               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 pe = cov_i - q_;
                 double alpha = (pe > 0.0) ? alphaPos : alphaNeg;
                 q_ += alpha * pe;
               } else {
                 pe = 0;
               }

               pes_[j] = pe;
               out_[j] = q_;
             }

             mark_run_complete();
           }
};

// 2D PE kernel: separate from DeltaKernel
struct Delta2Kernel : SequentialKernel {
  double qFast_ = NA_REAL;
  double qSlow_ = NA_REAL;
  double q_     = NA_REAL;

  // [trial][0 = fast PE, 1 = slow PE]
  std::vector<std::array<double, 2>> pes2_;

  Delta2Kernel() {}

  void run(const KernelParsView& kernel_pars,
           const Rcpp::NumericVector& covariate,
           const Rcpp::IntegerVector& comp_idx) override {

             int n_trials = comp_idx.size();

             out_.assign(n_trials, NA_REAL);
             pes2_.assign(n_trials, std::array<double,2>{NA_REAL, NA_REAL});

             const double* q0_col = kernel_pars.cols[0];
             const double* alphaFast_col = kernel_pars.cols[1];
             const double* propSlow_col = kernel_pars.cols[3];
             const double* dSwitch_col = kernel_pars.cols[4];

             int row0 = comp_idx[0];
             double q0 = q0_col[row0];
             qFast_ = qSlow_ = q_ = q0;

             for (int j = 0; j < n_trials; ++j) {
               int row = comp_idx[j];
               double cov_i     = covariate[row];
               double alphaFast = alphaFast_col[row];
               double propSlow  = propSlow_col[row];
               double dSwitch   = dSwitch_col[row];
               double alphaSlow = propSlow * alphaFast;

               double peFast = NA_REAL;
               double peSlow = NA_REAL;

               if (!Rcpp::NumericVector::is_na(cov_i)) {
                 peFast = cov_i - qFast_;
                 peSlow = cov_i - qSlow_;

                 qFast_ += alphaFast * peFast;
                 qSlow_ += alphaSlow * peSlow;

                 double diff = std::abs(qFast_ - qSlow_);
                 q_ = (diff > dSwitch) ? qFast_ : qSlow_;
               }

               pes2_[j][0] = peFast;
               pes2_[j][1] = peSlow;
               out_[j]     = q_;
             }

             mark_run_complete();
           }
};

// ---- Type mapping + factory ----

KernelType to_kernel_type(const Rcpp::String& k);

std::unique_ptr<BaseKernel> make_kernel(KernelType kt);





// Kernel handler
// class KernelHandler {
// public:
//   using Id = std::string;
//
//   KernelHandler() = default;
//
//   // Create a new kernel instance under a given ID (no run yet).
//   void create_instance(const Id& id, KernelType type) {
//     if (kernels_.count(id) > 0) {
//       Rcpp::stop("KernelHandler: ID '%s' already exists", id.c_str());
//     }
//     kernels_.emplace(id, make_kernel(type));
//   }
//
//   // Create and run in one step; return reference to output.
//   // If the ID already exists, throws.
//   const std::vector<double>& create_and_run(const Id& id,
//                                             KernelType type,
//                                             const Rcpp::NumericMatrix& kp,
//                                             const Rcpp::NumericVector& cov,
//                                             const Rcpp::IntegerVector& comp_idx) {
//     if (kernels_.count(id) > 0) {
//       Rcpp::stop("KernelHandler: ID '%s' already exists", id.c_str());
//     }
//     auto k = make_kernel(type);
//     k->reset();
//     k->run(kp, cov, comp_idx);
//
//     BaseKernel* raw = k.get();                 // keep raw pointer before move
//     kernels_.emplace(id, std::move(k));        // move into map
//     return raw->get_output();                  // safe: raw still valid
//   }
//
//   // Get cached output (no re-run).
//   const std::vector<double>& get_output(const Id& id) const {
//     const BaseKernel& k = get_kernel_const(id);
//     if (!k.has_run()) {
//       Rcpp::stop("KernelHandler: Kernel '%s' has not been run yet", id.c_str());
//     }
//     return k.get_output();
//   }
//
//   // Optional: get PEs for *1D* delta kernels only (SimpleDelta, Delta2LR)
//   Rcpp::NumericVector get_pes(const Id& id) const {
//     const BaseKernel& base = get_kernel_const(id);
//     if (!base.has_run()) {
//       Rcpp::stop("KernelHandler: Kernel '%s' has not been run yet", id.c_str());
//     }
//
//     const DeltaKernel* dk = dynamic_cast<const DeltaKernel*>(&base);
//     if (!dk) {
//       Rcpp::stop("KernelHandler: Kernel '%s' does not expose 1D PEs", id.c_str());
//     }
//
//     const std::vector<double>& pes_cpp = dk->get_pes();
//     Rcpp::NumericVector pes_R(pes_cpp.size());
//     for (int i = 0; i < (int)pes_cpp.size(); ++i)
//       pes_R[i] = pes_cpp[i];
//     return pes_R;
//   }
//
// private:
//   std::unordered_map<Id, std::unique_ptr<BaseKernel>> kernels_;
//
//   BaseKernel& get_kernel(const Id& id) {
//     auto it = kernels_.find(id);
//     if (it == kernels_.end()) {
//       Rcpp::stop("KernelHandler: No kernel with ID '%s'", id.c_str());
//     }
//     return *(it->second);
//   }
//
//   const BaseKernel& get_kernel_const(const Id& id) const {
//     auto it = kernels_.find(id);
//     if (it == kernels_.end()) {
//       Rcpp::stop("KernelHandler: No kernel with ID '%s'", id.c_str());
//     }
//     return *(it->second);
//   }
// };
