#include "TrendEngine.h"

// =============================================================================
// Rcpp boundary helpers — used only during TrendPlan construction
// =============================================================================

static std::string sexp_to_str(SEXP x) {
  if (Rf_isNull(x) || TYPEOF(x) != STRSXP || Rf_length(x) == 0) return "";
  return std::string(CHAR(STRING_ELT(x, 0)));
}

static std::vector<std::string> sexp_to_strvec(SEXP x) {
  std::vector<std::string> out;
  if (Rf_isNull(x) || TYPEOF(x) != STRSXP) return out;
  int n = Rf_length(x);
  out.reserve(n);
  for (int i = 0; i < n; ++i)
    out.emplace_back(CHAR(STRING_ELT(x, i)));
  return out;
}

static std::string list_str(const Rcpp::List& lst, const char* field) {
  if (!lst.containsElementNamed(field)) return "";
  return sexp_to_str(lst[field]);
}

static std::vector<std::string> list_strvec(const Rcpp::List& lst, const char* field) {
  if (!lst.containsElementNamed(field)) return {};
  return sexp_to_strvec(lst[field]);
}

static bool list_bool(const Rcpp::List& lst, const char* field, bool def = false) {
  if (!lst.containsElementNamed(field)) return def;
  SEXP x = lst[field];
  if (Rf_isNull(x)) return def;
  return Rcpp::as<bool>(x);
}

static int list_int(const Rcpp::List& lst, const char* field, int def = 0) {
  if (!lst.containsElementNamed(field)) return def;
  SEXP x = lst[field];
  if (Rf_isNull(x)) return def;
  return Rcpp::as<int>(x);
}

// =============================================================================
// KernelSpec construction helpers — plain C++ after data extraction
// =============================================================================

static void build_first_level(KernelSpec& ks, const Rcpp::DataFrame& data)
{
  const int n = data.nrows();
  if (n <= 0) Rf_error("build_first_level: data has zero rows");

  ks.first_level.assign(n, true);

  if (ks.has_at) {
    if (!data.containsElementNamed(ks.at.c_str()))
      Rf_error("build_first_level: data has no column '%s'", ks.at.c_str());
    SEXP at_col = data[ks.at.c_str()];
    if (!Rf_inherits(at_col, "factor"))
      Rf_error("'at' column '%s' must be a factor", ks.at.c_str());
    if (Rf_length(at_col) != n)
      Rf_error("'at' column '%s' has wrong length", ks.at.c_str());
    const int* f = INTEGER(at_col);
    for (int i = 0; i < n; ++i)
      ks.first_level[i] = (f[i] == 1);
  }

  ks.expand_idx.assign(n, 0);
  int count = 0;
  for (int i = 0; i < n; ++i) {
    if (ks.first_level[i]) ++count;
    ks.expand_idx[i] = count;
  }
  if (count == 0)
    Rf_error("build_first_level: no rows with first 'at' level found");
  for (int i = 0; i < n; ++i)
    if (ks.expand_idx[i] == 0)
      Rf_error("build_first_level: rows before first 'at' level");

  ks.comp_index.clear();
  ks.comp_index.reserve(count);
  for (int i = 0; i < n; ++i)
    if (ks.first_level[i]) ks.comp_index.push_back(i);
}

static void build_kernel_input(KernelSpec& ks, const Rcpp::DataFrame& data)
{
  const int n_cov = (int)ks.cov_names.size();
  const int n_par = (int)ks.par_input.size();
  const int n_col = n_cov + n_par;
  const int n_row = data.nrows();

  if (n_col == 0)
    Rf_error("KernelSpec '%s': needs at least one cov_name or par_input",
             ks.kernel_id.c_str());

  ks.kernel_input = Mat(n_row, n_col);
  ks.covariate_indices.clear();
  ks.par_input_indices.clear();

  for (int i = 0; i < n_cov; ++i) {
    const char* cn = ks.cov_names[i].c_str();
    if (!data.containsElementNamed(cn))
      Rf_error("KernelSpec '%s': data has no column '%s'",
               ks.kernel_id.c_str(), cn);
    SEXP col = data[cn];
    if (TYPEOF(col) != REALSXP && TYPEOF(col) != INTSXP)
      Rf_error("KernelSpec '%s': covariate column '%s' must be numeric or integer",
               ks.kernel_id.c_str(), cn);

    double* dst = ks.kernel_input.colptr(i);
    if (TYPEOF(col) == REALSXP) {
      const double* src = REAL(col);
      std::copy(src, src + n_row, dst);
    } else {
      const int* src = INTEGER(col);
      for (int j = 0; j < n_row; ++j)
        dst[j] = (src[j] == NA_INTEGER) ? NA_REAL : static_cast<double>(src[j]);
    }
    ks.covariate_indices.push_back(i);
  }

  for (int i = 0; i < n_par; ++i)
    ks.par_input_indices.push_back(n_cov + i);
  // par_input columns remain zero-initialised; filled per-particle at runtime
}

static void build_kernel_args(KernelSpec& ks,
                              const Rcpp::List& k_lst,
                              const Rcpp::DataFrame& data)
{
  ks.q_reset_col.clear();

  if (!k_lst.containsElementNamed("kernel_args")) { ks.build_kernel_args(); return; }
  SEXP ka_sexp = k_lst["kernel_args"];
  if (Rf_isNull(ka_sexp))                          { ks.build_kernel_args(); return; }

  Rcpp::List ka(ka_sexp);
  if (ka.containsElementNamed("q_reset_column")) {
    SEXP cns = ka["q_reset_column"];
    if (!Rf_isNull(cns)) {
      std::string col_name = sexp_to_str(cns);
      if (!data.containsElementNamed(col_name.c_str()))
        Rf_error("kernel_args$q_reset_column: column '%s' not found", col_name.c_str());
      SEXP col = data[col_name.c_str()];
      int n = data.nrows();
      ks.q_reset_col.resize(n);
      if      (TYPEOF(col) == LGLSXP) { const int* p = LOGICAL(col); std::copy(p, p+n, ks.q_reset_col.data()); }
      else if (TYPEOF(col) == INTSXP)  { const int* p = INTEGER(col); std::copy(p, p+n, ks.q_reset_col.data()); }
      else Rf_error("kernel_args$q_reset_column: column '%s' must be logical or integer",
                    col_name.c_str());
      if ((int)ks.q_reset_col.size() != n)
        Rf_error("kernel_args$q_reset_column: wrong length");
    }
  }
  if (ka.containsElementNamed("grid_res")) {
    SEXP gr = ka["grid_res"];
    if (!Rf_isNull(gr)) ks.kernel_args.grid_res = Rcpp::as<int>(gr);
  }
  ks.build_kernel_args();
}

static void build_covariate_coding(BaseSpec& bs,
                                   const Rcpp::List& b_lst,
                                   const Rcpp::DataFrame& data,
                                   const Rcpp::List& data_covcoding)
{
  bs.has_covariate_coding = false;
  bs.covariate_coding.clear();

  if (!b_lst.containsElementNamed("coding") || Rf_isNull(b_lst["coding"])) return;
  if (data_covcoding.size() == 0)
    Rf_error("BaseSpec '%s': 'coding' specified but data has no 'covariate_coding' attribute",
             bs.target_parameter.c_str());

  Rcpp::List coding_spec(b_lst["coding"]);
  Rcpp::CharacterVector coding_names = coding_spec.names();
  if (coding_names.size() != coding_spec.size())
    Rf_error("BaseSpec '%s': 'coding' list must be named", bs.target_parameter.c_str());
  if (coding_spec.size() != 1)
    Rf_error("BaseSpec '%s': 'coding' must contain exactly one entry",
             bs.target_parameter.c_str());

  Rcpp::CharacterVector data_coding_names = data_covcoding.names();
  const int T = data.nrows();

  std::string coding_nm = Rcpp::as<std::string>(coding_names[0]);

  int idx = -1;
  for (int j = 0; j < data_coding_names.size(); ++j)
    if (Rcpp::as<std::string>(data_coding_names[j]) == coding_nm) { idx = j; break; }
  if (idx < 0)
    Rf_error("BaseSpec '%s': covariate_coding '%s' not found in data",
             bs.target_parameter.c_str(), coding_nm.c_str());

  SEXP mat_sexp = data_covcoding[idx];
  if (Rf_ncols(mat_sexp) == 0)
    Rf_error("BaseSpec '%s': covariate_coding['%s'] has zero columns",
             bs.target_parameter.c_str(), coding_nm.c_str());
  if (Rf_nrows(mat_sexp) != T)
    Rf_error("BaseSpec '%s': covariate_coding['%s'] has wrong nrow",
             bs.target_parameter.c_str(), coding_nm.c_str());

  bs.covariate_coding.push_back(Mat::from_sexp(mat_sexp));
  bs.has_covariate_coding = true;
}

// =============================================================================
// TrendPlan constructor — Rcpp boundary
// =============================================================================

TrendPlan::TrendPlan(const Rcpp::List& trend, const Rcpp::DataFrame& data)
{
  // covariate_coding attribute on data
  Rcpp::List data_covcoding;
  if (data.hasAttribute("covariate_coding")) {
    SEXP cm = data.attr("covariate_coding");
    if (!Rf_isNull(cm)) data_covcoding = Rcpp::List(cm);
  }

  // ------------------------------------------------------------------
  // 1. Parse kernels
  // ------------------------------------------------------------------
  if (!trend.containsElementNamed("kernels"))
    Rf_error("TrendPlan: trend list must have a 'kernels' element");

  Rcpp::List k_list = Rcpp::as<Rcpp::List>(trend["kernels"]);
  Rcpp::CharacterVector k_names = k_list.names();

  for (int i = 0; i < k_list.size(); ++i) {
    Rcpp::List k_lst = Rcpp::as<Rcpp::List>(k_list[i]);
    KernelSpec ks;

    ks.kernel_id   = Rcpp::as<std::string>(k_names[i]);
    ks.kernel      = list_str(k_lst, "type");
    ks.kernel_type = to_kernel_type(Rcpp::String(ks.kernel));
    ks.phase       = parse_phase(list_str(k_lst, "phase"));
    ks.sequential  = list_bool(k_lst, "sequential", false);
    ks.cov_names   = list_strvec(k_lst, "cov_names");
    ks.par_input   = list_strvec(k_lst, "par_input");
    ks.pnames      = list_strvec(k_lst, "pnames");

    if (k_lst.containsElementNamed("at") && !Rf_isNull(k_lst["at"])) {
      ks.has_at = true;
      ks.at     = list_str(k_lst, "at");
    }

    if (ks.kernel_type == KernelType::Custom) {
      if (!k_lst.containsElementNamed("kernel_pointer"))
        Rf_error("Custom kernel '%s' requires 'kernel_pointer' field",
                 ks.kernel_id.c_str());
      SEXP ptr = k_lst["kernel_pointer"];
      if (Rf_isNull(ptr))
        Rf_error("Custom kernel '%s': 'kernel_pointer' is NULL", ks.kernel_id.c_str());
      ks.custom_fun = ptr;
    }

    build_kernel_args(ks, k_lst, data);
    build_kernel_input(ks, data);
    build_first_level(ks, data);

    // populate param sets
    for (const auto& pn : ks.pnames) {
      all_trend_params.insert(pn);
      switch (ks.phase) {
      case TrendPhase::Premap:        premap_params.insert(pn);        break;
      case TrendPhase::Pretransform:  pretransform_params.insert(pn);  break;
      case TrendPhase::Posttransform: posttransform_params.insert(pn); break;
      }
    }

    kernels.emplace(ks.kernel_id, std::move(ks));
  }

  // ------------------------------------------------------------------
  // 2. Parse bases
  // ------------------------------------------------------------------
  if (!trend.containsElementNamed("bases"))
    Rf_error("TrendPlan: trend list must have a 'bases' element");

  Rcpp::List b_list = Rcpp::as<Rcpp::List>(trend["bases"]);

  for (int i = 0; i < b_list.size(); ++i) {
    Rcpp::List b_lst = Rcpp::as<Rcpp::List>(b_list[i]);
    BaseSpec bs;

    bs.base_type        = list_str(b_lst, "type");
    bs.target_parameter = list_str(b_lst, "target_parameter");
    bs.kernel_id        = list_str(b_lst, "kernel_id");
    bs.kernel_output    = list_int(b_lst, "kernel_output", 1);
    bs.phase            = parse_phase(list_str(b_lst, "phase"));
    bs.pnames           = list_strvec(b_lst, "pnames");

    if (kernels.find(bs.kernel_id) == kernels.end())
      Rf_error("BaseSpec '%s': kernel_id '%s' not found",
               bs.target_parameter.c_str(), bs.kernel_id.c_str());

    build_covariate_coding(bs, b_lst, data, data_covcoding);

    for (const auto& pn : bs.pnames) {
      all_trend_params.insert(pn);
      switch (bs.phase) {
      case TrendPhase::Premap:        premap_params.insert(pn);        break;
      case TrendPhase::Pretransform:  pretransform_params.insert(pn);  break;
      case TrendPhase::Posttransform: posttransform_params.insert(pn); break;
      }
    }

    switch (bs.phase) {
    case TrendPhase::Premap:        premap_bases.push_back(std::move(bs));        break;
    case TrendPhase::Pretransform:  pretransform_bases.push_back(std::move(bs));  break;
    case TrendPhase::Posttransform: posttransform_bases.push_back(std::move(bs)); break;
    }
  }
}

// =============================================================================
// TrendPlan::premap_design_mask  — Rcpp boundary
// =============================================================================

Rcpp::LogicalVector TrendPlan::premap_design_mask(const Rcpp::List& designs) const
{
  const int n = designs.size();
  Rcpp::LogicalVector mask(n, false);
  Rcpp::CharacterVector dnames = designs.names();
  for (int i = 0; i < n; ++i) {
    std::string nm = Rcpp::as<std::string>(dnames[i]);
    mask[i] = (premap_params.count(nm) > 0);
  }
  return mask;
}

// =============================================================================
// TrendRuntime constructor
// =============================================================================

TrendRuntime::TrendRuntime(const TrendPlan& plan_) : plan(&plan_)
{
  for (const auto& kv : plan->kernels) {
    const KernelSpec& ks = kv.second;
    KernelRuntime k_rt;
    k_rt.spec       = &ks;

    // determine number of kernels to create based on kernel's arity and number of covariates
    KernelMeta meta = kernel_meta(ks.kernel_type);
    const bool variadic = (meta.input_arity == -1 ||
                           ks.kernel_type == KernelType::Custom);

    int n_slots = variadic ? 1                             // variadic or custom: one kernel sees all columns
    : (int)(ks.cov_names.size() + ks.par_input.size());    // arity-1: one kernel per covariate

    for (int s = 0; s < n_slots; ++s) {
      auto kptr = make_kernel(ks.kernel_type, ks.custom_fun);
      kptr->set_kernel_args(ks.kernel_args);
      k_rt.kernel_ptrs.push_back(std::move(kptr));
    }

    // fill data covariates, pre-allocate par_input memory
    if (!variadic) {
      const int n = ks.kernel_input.nrow;
      k_rt.slot_inputs.reserve(n_slots);

      // covariate slots: copy data in once, never touched again
      for (int s = 0; s < (int)ks.covariate_indices.size(); ++s) {
        Mat buf(n, 1);
        std::copy(ks.kernel_input.colptr(ks.covariate_indices[s]),
                  ks.kernel_input.colptr(ks.covariate_indices[s]) + n,
                  buf.colptr(0));
        k_rt.slot_inputs.push_back(std::move(buf));
      }

      // par_input slots: allocate zero buffer, record which par_input feeds it
      for (int i = 0; i < (int)ks.par_input_indices.size(); ++i) {
        int slot_idx = (int)ks.covariate_indices.size() + i;
        k_rt.slot_inputs.push_back(Mat(n, 1));  // zero-init, filled per particle
        k_rt.par_input_slot_indices.push_back(slot_idx);
        k_rt.par_input_param_col.push_back(i);  // index into ks.par_input
      }
    }
    kernels.emplace(ks.kernel_id, std::move(k_rt));
  }

  auto build_base_rts = [&](const std::vector<BaseSpec>& specs,
                            std::vector<BaseRuntime>& out) {
    out.reserve(specs.size());
    for (const auto& bs : specs) {
      BaseRuntime b_rt;
      b_rt.spec      = &bs;
      b_rt.kernel_rt = &kernels.at(bs.kernel_id);
      out.push_back(std::move(b_rt));
    }
  };

  build_base_rts(plan->premap_bases,        premap_bases);
  build_base_rts(plan->pretransform_bases,  pretransform_bases);
  build_base_rts(plan->posttransform_bases, posttransform_bases);
}

// =============================================================================
// TrendRuntime::bind_all_to_paramtable
// =============================================================================

void TrendRuntime::bind_all_to_paramtable(const ParamTable& pt)
{
  for (auto& kv : kernels) {
    KernelRuntime& k_rt = kv.second;
    const KernelSpec& ks = *k_rt.spec;
    k_rt.kernel_par_indices.clear();
    k_rt.kernel_par_indices.reserve(ks.pnames.size());
    for (const auto& pn : ks.pnames)
      k_rt.kernel_par_indices.push_back(pt.base_index_for(pn));
    // validate par_input names
    for (const auto& pn : ks.par_input)
      (void)pt.base_index_for(pn);
  }

  auto bind_bases = [&](std::vector<BaseRuntime>& base_rts) {
    for (auto& b_rt : base_rts) {
      const BaseSpec& bs = *b_rt.spec;
      b_rt.base_par_indices.clear();
      b_rt.base_par_indices.reserve(bs.pnames.size());
      for (const auto& pn : bs.pnames)
        b_rt.base_par_indices.push_back(pt.base_index_for(pn));
      b_rt.kernel_rt = &kernels.at(bs.kernel_id);
    }
  };

  bind_bases(premap_bases);
  bind_bases(pretransform_bases);
  bind_bases(posttransform_bases);
}

// =============================================================================
// TrendRuntime::reset_all_kernels
// =============================================================================

void TrendRuntime::reset_all_kernels() {
  for (auto& kv : kernels)
    for (auto& kptr : kv.second.kernel_ptrs)
      kptr->reset();
}

// =============================================================================
// TrendRuntime::run_kernel  (private)
// =============================================================================

void TrendRuntime::run_kernel(KernelRuntime& k_rt, ParamTable& pt)
{
  const KernelSpec& ks   = *k_rt.spec;
  const int         n    = pt.n_trials;
  const bool variadic = k_rt.is_variadic();

  if (variadic) {
    // fill par_input columns into shared matrix in-place
    for (int i = 0; i < (int)ks.par_input_indices.size(); ++i) {
      const double* src = pt.column_by_name_ptr(ks.par_input[i]);
      double*       dst = ks.kernel_input.colptr(ks.par_input_indices[i]);
      std::copy(src, src + n, dst);
    }
    auto& kptr = k_rt.kernel_ptrs[0];
    kptr->reset();
    kptr->run(make_kernel_pars_view(pt, k_rt.kernel_par_indices),
              ks.kernel_input, ks.comp_index);
    if (ks.has_at) {
      kptr->set_expand_idx(ks.expand_idx);
      kptr->do_expand(ks.expand_idx);
    }
  } else {
    // overwrite par_input slot buffers in-place — no allocation
    for (int j = 0; j < (int)k_rt.par_input_slot_indices.size(); ++j) {
      int slot_idx  = k_rt.par_input_slot_indices[j];
      int par_i     = k_rt.par_input_param_col[j];
      const double* src = pt.column_by_name_ptr(ks.par_input[par_i]);
      double*       dst = k_rt.slot_inputs[slot_idx].colptr(0);
      std::copy(src, src + n, dst);
    }

    // run each slot kernel against its pre-allocated buffer
    for (int s = 0; s < k_rt.n_slots(); ++s) {
      auto& kptr = k_rt.kernel_ptrs[s];
      kptr->reset();
      kptr->run(make_kernel_pars_view(pt, k_rt.kernel_par_indices),
                k_rt.slot_inputs[s], ks.comp_index);
      if (ks.has_at) {
        kptr->set_expand_idx(ks.expand_idx);
        kptr->do_expand(ks.expand_idx);
      }
    }
  }
}

// =============================================================================
// TrendRuntime::apply_base
// =============================================================================

void TrendRuntime::apply_base(BaseRuntime& base_rt, ParamTable& pt)
{
  KernelRuntime&  k_rt = *base_rt.kernel_rt;
  const BaseSpec& bs   = *base_rt.spec;

  if (!k_rt.kernel_ptrs[0]->has_run())
    run_kernel(k_rt, pt);

  const int n = pt.n_trials;

  const int target_idx = pt.base_index_for(bs.target_parameter);
  double*   target_col = &pt.base(0, target_idx);

  const bool base_is_lin      = (bs.base_type == "lin");
  const bool base_is_centered = (bs.base_type == "centered");
  const bool base_is_identity = (bs.base_type == "identity");
  const bool needs_base_par   = base_is_lin || base_is_centered;
  const double center_offset = base_is_centered ? 0.5 : 0.0;

  // single base par pointer — same par applied to all slots
  const double* base_ptr = needs_base_par
  ? &pt.base(0, base_rt.base_par_indices[0])
    : nullptr;

  const bool variadic = k_rt.is_variadic();
  const int  n_slots  = k_rt.n_slots();

  const bool has_coding = bs.has_covariate_coding && !bs.covariate_coding.empty();
  const Mat* coding     = has_coding ? &bs.covariate_coding[0] : nullptr;

  // -----------------------------------------------------------------------
  // Variadic kernel: single kernel, output is [n x n_cols]
  // -----------------------------------------------------------------------
  if (variadic) {
    KernelOutput ko = k_rt.kernel_ptrs[0]->get_output_stream(bs.kernel_output);

    if (has_coding && needs_base_par) {
      for (int col = 0; col < ko.n_cols; ++col) {
        const double* q_col = ko.data + col * ko.n_rows;
        const double* m_col = coding->colptr(col < coding->ncol ? col : 0);
#pragma omp simd
        for (int r = 0; r < n; ++r) {
          double q = q_col[r] - center_offset;
          target_col[r] += q * base_ptr[r] * m_col[r];
        }
      }
    } else if (needs_base_par) {
      for (int col = 0; col < ko.n_cols; ++col) {
        const double* q_col = ko.data + col * ko.n_rows;
#pragma omp simd
        for (int r = 0; r < n; ++r) {
          double q = q_col[r] - center_offset;
          target_col[r] += q * base_ptr[r];
        }
      }
    } else if (base_is_identity) {
      // identity replaces rather than accumulates — intentional
      for (int col = 0; col < ko.n_cols; ++col) {
        const double* q_col = ko.data + col * ko.n_rows;
#pragma omp simd
        for (int r = 0; r < n; ++r)
          target_col[r] = q_col[r];
      }
    } else {
      for (int col = 0; col < ko.n_cols; ++col) {
        const double* q_col = ko.data + col * ko.n_rows;
#pragma omp simd
        for (int r = 0; r < n; ++r)
          target_col[r] += q_col[r];
      }
    }

    // -----------------------------------------------------------------------
    // Arity-1: one kernel per slot, map column s scales slot s
    // -----------------------------------------------------------------------
  } else {
    for (int s = 0; s < n_slots; ++s) {
      KernelOutput  ko    = k_rt.kernel_ptrs[s]->get_output_stream(bs.kernel_output);
      const double* q_col = ko.data;

      if (has_coding && needs_base_par) {
        const double* m_col = coding->colptr(s);
#pragma omp simd
        for (int r = 0; r < n; ++r) {
          double q = q_col[r] - center_offset;
          target_col[r] += q * base_ptr[r] * m_col[r];
        }
      } else if (needs_base_par) {
#pragma omp simd
        for (int r = 0; r < n; ++r) {
          double q = q_col[r] - center_offset;
          target_col[r] += q * base_ptr[r];
        }
      } else if (base_is_identity) {
        // identity replaces rather than accumulates — intentional
#pragma omp simd
        for (int r = 0; r < n; ++r)
          target_col[r] = q_col[r];
      } else {
#pragma omp simd
        for (int r = 0; r < n; ++r)
          target_col[r] += q_col[r];
      }
    }
  }
}



// =============================================================================
// TrendRuntime::all_kernel_outputs  — Rcpp boundary (diagnostic only)
// =============================================================================

Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt)
{
  return all_kernel_outputs(pt, std::vector<int>{1});
}

Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt,
                                                     const std::vector<int>& codes)
{
  const int n = pt.n_trials;

  // ensure all kernels have run
  for (auto& kv : kernels)
    if (!kv.second.kernel_ptrs[0]->has_run())
      run_kernel(kv.second, pt);

  // count total output columns across all slots
  int n_cols_total = 0;
  for (const auto& kv : kernels) {
    const KernelRuntime& k_rt = kv.second;
    for (int s = 0; s < k_rt.n_slots(); ++s)
      for (int code : codes) {
        if (k_rt.kernel_ptrs[s]->has_output_stream(code))
          n_cols_total += k_rt.kernel_ptrs[s]->get_output_stream(code).n_cols;
      }
  }

  Rcpp::NumericMatrix out(n, n_cols_total);
  Rcpp::CharacterVector cn(n_cols_total);
  int col = 0;

  for (auto& kv : kernels) {
    KernelRuntime&    k_rt = kv.second;
    const KernelSpec& ks   = *k_rt.spec;

    // unified input names in slot order: covariates first, then par_inputs
    std::vector<std::string> input_names;
    input_names.insert(input_names.end(), ks.cov_names.begin(),  ks.cov_names.end());
    input_names.insert(input_names.end(), ks.par_input.begin(),  ks.par_input.end());

    for (int s = 0; s < k_rt.n_slots(); ++s) {
      for (int code : codes) {
        if (!k_rt.kernel_ptrs[s]->has_output_stream(code)) continue;
        KernelOutput ko     = k_rt.kernel_ptrs[s]->get_output_stream(code);
        std::string  suffix = k_rt.kernel_ptrs[s]->output_stream_name(code);

        if (ko.n_rows != n)
          Rf_error("all_kernel_outputs '%s': stream rows != n_trials",
                   ks.kernel_id.c_str());

        for (int c = 0; c < ko.n_cols; ++c) {
          const double* src = ko.data + c * ko.n_rows;
          std::copy(src, src + n, &out(0, col));

          // arity-1: slot s maps directly to input_names[s]
          // variadic: column c maps to input_names[c] (all inputs fed as one matrix)
          int name_idx      = k_rt.is_variadic() ? c : s;
          std::string iname = (name_idx < (int)input_names.size())
            ? input_names[name_idx]
          : std::to_string(name_idx + 1);
          cn[col] = ks.kernel_id + "." + iname + "." + suffix;
          ++col;
        }
      }
    }
  }

  Rcpp::colnames(out) = cn;
  return out;
}

