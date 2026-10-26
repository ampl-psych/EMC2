#include <Rcpp.h>
#include <algorithm>
#include "CensorSpec.h"
#include "composite_functions.h"

using namespace Rcpp;

// One direct time integral for all independent race models. Ordinary
// probabilities use the existing cubature; log scaling handles underflow.
double CensorSpec::log_interval(int base,int winner,double lower,double upper,
                               const std::vector<int>& is_ok) const {
  if (!(upper>lower)) return neg_inf();
  if (setup->log_response_interval)
    return setup->log_response_interval(*pt,setup->spec,base,winner+1,lower,upper);
  const int win=base+winner;
  if (!participating[win]) return neg_inf();
  lower=std::max(lower,pt->base(win,setup->spec.col_t0));
  if (!(upper>lower)) return neg_inf();
  winner_row[0]=win;
  active_rows.clear();
  breaks.clear();
  breaks.push_back(lower);
  for (int k=0;k<n_acc;++k) {
    const int row=base+k;
    if (!participating[row]) continue;
    if (row!=win) active_rows.push_back(row);
    const double onset=pt->base(row,setup->spec.col_t0);
    if (onset>lower && onset<upper) breaks.push_back(onset);
  }
  breaks.push_back(upper);
  std::sort(breaks.begin(),breaks.end());
  breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
  auto density=[&](double t,bool log_p) {
    bound_work[win]=t;
    for (int row:active_rows) bound_work[row]=t;
    setup->fill_both(bound_work.data(),*pt,setup->spec,winner_row,active_rows,
                     log_p_lower.data(),*scratch,is_ok,log_p);
    double value=log_p_lower[win];
    for (int row:active_rows)
      if (log_p) value+=log_p_lower[row]; else value*=log_p_lower[row];
    return value;
  };
  double result=neg_inf();
  for (std::size_t segment=1;segment<breaks.size();++segment) {
    const double lo=breaks[segment-1],hi=breaks[segment];
    const bool infinite=is_pos_inf(hi);
    const double scale=infinite ? std::max(1.0,std::fabs(lo)) : hi-lo;
    auto probability=[&](double u) {
      const double denominator=infinite ? 1-u : 1.0;
      if (denominator==0.0) return 0.0;
      return density(lo+scale*u/denominator,false)*scale/(denominator*denominator);
    };
    using Function=decltype(probability);
    auto integrand=[](unsigned,const double* x,void* p,unsigned,double* out) -> int {
      out[0]=(*static_cast<Function*>(p))(x[0]);return 0;
    };
    double a=0,b=1,value=0,error=0;
    constexpr double tolerance=1e-8;
    hcubature(integrand,&probability,1,&a,&b,6000,0,tolerance,&value,&error);
    double lp;
    if (is_finite(value) && value>1e-280 && error<=tolerance*value) lp=std::log(value);
    else {
      auto log_probability=[&](double u) {
        if (infinite && u==1.0) return neg_inf();
        const double t=infinite ? lo+scale*u/(1-u) : lo+scale*u;
        return density(t,true)+std::log(scale)-(infinite ? 2*std::log1p(-u) : 0.0);
      };
      lp=log_integral_unit(log_probability,-1.0,tolerance);
    }
    result=log_sum_exp(result,lp);
  }
  return result;
}

// Censored likelihoods remain on the log scale through truncation normalization.
void CensorSpec::fill_censored_rows(const TruncSpec& trunc,
                                    const std::vector<int>& is_ok,
                                    std::vector<double>& ll_trial) const
{
  if (setup->log_response_interval) {
    auto fill=[&](const std::vector<int>& idx,bool left,bool right) {
      for (int trial:idx) {
        const int base=trial*n_acc;
        if (!is_ok[base]) continue;
        const double lt=trunc.LT.empty() ? 0.0 : trunc.LT[base];
        const double ut=trunc.UT.empty() ? pos_inf() : trunc.UT[base];
        double lp=neg_inf();
        if (left) lp=log_interval(base,-1,lt,std::min(LC[base],ut),is_ok);
        if (right) lp=log_sum_exp(lp,log_interval(base,-1,std::max(UC[base],lt),ut,is_ok));
        ll_trial[trial]=lp;
      }
    };
    fill(idx_L_tr,true,false);
    fill(idx_U_tr,false,true);
    fill(idx_B_tr,true,true);
  } else {
    std::fill(log_p_lower.begin(), log_p_lower.end(), 0.0);
    std::fill(log_p_upper.begin(), log_p_upper.end(), 0.0);
    setup->fill_survivor(idx_L, LC, *pt, setup->spec, log_p_lower.data(), *scratch, is_ok, true);
    setup->fill_survivor(idx_U, UC, *pt, setup->spec, log_p_upper.data(), *scratch, is_ok, true);

    auto race_survival = [&](const std::vector<double>& p,int trial) {
      double sum=0.0;
      for (int k=0; k<n_acc; ++k) {
        const int row=trial*n_acc+k;
        if (participating[row]) sum+=p[row];
      }
      return sum;
    };
    for (int row : idx_L_tr) {
      if (!is_ok[row*n_acc]) continue;
      ll_trial[row]=log_diff_exp(trunc.log_S_lower[row],race_survival(log_p_lower,row));
    }
    for (int row : idx_U_tr) {
      if (!is_ok[row*n_acc]) continue;
      ll_trial[row]=log_diff_exp(race_survival(log_p_upper,row),trunc.log_S_upper[row]);
    }
    for (int row : idx_B_tr) {
      if (!is_ok[row*n_acc]) continue;
      ll_trial[row]=log_sum_exp(
        log_diff_exp(trunc.log_S_lower[row],race_survival(log_p_lower,row)),
        log_diff_exp(race_survival(log_p_upper,row),trunc.log_S_upper[row]));
    }
  }

  // Known responses and silent accumulators share the same interval evaluator.
  auto fill_known = [&](const std::vector<int>& idx, const std::vector<int>& winner,
                        bool lower, bool upper) {
    for (std::size_t j=0; j<idx.size(); ++j) {
      const int base=idx[j], trial=base/n_acc;
      if (!is_ok[base]) continue;
      const double lt=trunc.LT.empty() ? 0.0 : trunc.LT[base];
      const double ut=trunc.UT.empty() ? pos_inf() : trunc.UT[base];
      double lp=neg_inf();
      if (lower) lp=log_interval(base,winner[j],lt,std::min(LC[base],ut),is_ok);
      if (upper) lp=log_sum_exp(lp,log_interval(base,winner[j],std::max(UC[base],lt),ut,is_ok));
      ll_trial[trial]=lp;
    }
  };
  fill_known(idx_L_known,winner_L_known,true,false);
  fill_known(idx_U_known,winner_U_known,false,true);
  fill_known(idx_B_known,winner_B_known,true,true);
  for (std::size_t j=0; j<idx_silent.size(); ++j) {
    const int base=idx_silent[j], trial=base/n_acc;
    if (!is_ok[base]) continue;
    const double lt=trunc.LT.empty() ? 0.0 : trunc.LT[base];
    const double ut=trunc.UT.empty() ? pos_inf() : trunc.UT[base];
    ll_trial[trial]=log_sum_exp(ll_trial[trial],
      log_interval(base,winner_silent[j],lt,std::min(UC[base+winner_silent[j]],ut),is_ok));
  }
}

// -----------------------------------------------------------------------------
// make_censor_spec()
//
// missingness codes:
//   1 = lower-censored
//   2 = upper-censored
//   3 = both-censored
//
// For each censored trial, data$R determines the path:
//   NA  → analytical (no known response)
//   1,2,... → integration (known response; winner = R - 1, 0-based)
// -----------------------------------------------------------------------------
CensorSpec make_censor_spec(const DataFrame& data,
                            int n_trials,
                            int n_acc,
                            const std::vector<bool>& participating,
                            const RaceModelSetup* setup,
                            const ParamTable& pt,
                            RaceScratch& scratch)
{
  CensorSpec censor;
  censor.n_trials = n_trials;
  censor.n_acc    = n_acc;
  censor.n_rows   = n_trials * n_acc;
  censor.setup    = setup;
  censor.participating = participating;
  censor.pt       = &pt;
  censor.scratch  = &scratch;

  std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(data.names());
  const bool has_LC_col      = std::find(names.begin(), names.end(), "LC")          != names.end();
  const bool has_UC_col      = std::find(names.begin(), names.end(), "UC")          != names.end();
  const bool has_missingness = std::find(names.begin(), names.end(), "missingness") != names.end();

  if (!has_missingness) return censor;

  IntegerVector missingness_tmp = data["missingness"];
  const int*    missingness_ptr = INTEGER(missingness_tmp);

  IntegerVector R_tmp = data["R"];
  const int*    R_ptr = INTEGER(R_tmp);

  NumericVector lc_tmp, uc_tmp;
  if (has_LC_col) {
    lc_tmp = data["LC"];
    censor.LC.assign(REAL(lc_tmp), REAL(lc_tmp) + censor.n_rows);
    censor.idx_L.reserve(censor.n_rows);
    censor.idx_L_tr.reserve(censor.n_trials);
    censor.idx_L_known.reserve(censor.n_trials);
    censor.winner_L_known.reserve(censor.n_trials);
  }
  if (has_UC_col) {
    uc_tmp = data["UC"];
    censor.UC.assign(REAL(uc_tmp), REAL(uc_tmp) + censor.n_rows);
    censor.idx_U.reserve(censor.n_rows);
    censor.idx_U_tr.reserve(censor.n_trials);
    censor.idx_B_tr.reserve(censor.n_trials);
    censor.idx_U_known.reserve(censor.n_trials);
    censor.idx_B_known.reserve(censor.n_trials);
    censor.winner_U_known.reserve(censor.n_trials);
    censor.winner_B_known.reserve(censor.n_trials);
  }

  // Row-level loop: populate per-accumulator index lists (analytical cases only)
  for (int i = 0; i < censor.n_rows; ++i) {
    switch (missingness_ptr[i]) {
    case 1: censor.idx_L.push_back(i); break;
    case 2: censor.idx_U.push_back(i); break;
    case 3: censor.idx_L.push_back(i); censor.idx_U.push_back(i); break;
    case 4: censor.idx_U.push_back(i); break; // silent accumulator contributes to S_race(UC)
    default: break;
    }
  }

  censor.idx_silent.reserve(n_trials);
  censor.winner_silent.reserve(n_trials);

  // Trial-level loop (step by n_acc): populate trial index lists
  for (int base = 0; base < censor.n_rows; base += n_acc) {
    // Aggregate missingness across all accumulator rows for this trial
    bool has_L = false, has_U = false;
    bool has_silent = false;
    int  R_val = NA_INTEGER;
    for (int k = 0; k < n_acc; ++k) {
      const int m = missingness_ptr[base + k];
      if (m == 1 || m == 3) has_L = true;
      if (m == 2 || m == 3) has_U = true;
      if (m == 4)           has_silent   = true;
      // Take first non-NA R value as the winner for regular censoring with known response
      // (cases without silent accumulators)
      if (R_val == NA_INTEGER && R_ptr[base + k] != NA_INTEGER)
        R_val = R_ptr[base + k];
    }

    // For cases with silent accumulators, we need to mark the accumulators with missingness == 4 as the winner
    if (has_silent) {
      // Silent accumulators contribute S_race(UC) via idx_U_tr
      censor.idx_U_tr.push_back(base / n_acc);
      for (int k = 0; k < n_acc; ++k) {
        if (missingness_ptr[base + k] == 4) {
          censor.idx_silent.push_back(base);
          censor.winner_silent.push_back(k);
          }
      }
      continue;  // skip regular miss aggregation for this trial
    }

    const int miss = (has_L ? 1 : 0) + (has_U ? 2 : 0);
    if (miss == 0) continue;

    const bool known = (R_val != NA_INTEGER);
    const int winner = known ? (R_val - 1) : -1;

    switch (miss) {
    case 1:
      if (known) { censor.idx_L_known.push_back(base); censor.winner_L_known.push_back(winner); }
      else         censor.idx_L_tr.push_back(base / n_acc);
      break;
    case 2:
      if (known) { censor.idx_U_known.push_back(base); censor.winner_U_known.push_back(winner); }
      else         censor.idx_U_tr.push_back(base / n_acc);
      break;
    case 3:
      if (known) { censor.idx_B_known.push_back(base); censor.winner_B_known.push_back(winner); }
      else         censor.idx_B_tr.push_back(base / n_acc);
      break;
    default: break;
    }
  }


  // Validate consistency
  if (!censor.idx_L.empty() && censor.LC.empty())
    Rcpp::stop("missingness contains lower-censored rows (1 or 3) but no LC column found in data");
  if (!censor.idx_U.empty() && censor.UC.empty())
    Rcpp::stop("missingness contains upper-censored rows (2, 3 or 4) but no UC column found in data");

  std::vector<bool> unknown(n_trials,false);
  for (int trial:censor.idx_L_tr) unknown[trial]=true;
  for (int trial:censor.idx_U_tr) unknown[trial]=true;
  for (int trial:censor.idx_B_tr) unknown[trial]=true;
  auto known=[&](int row) { return !unknown[row/n_acc]; };
  censor.idx_L.erase(std::remove_if(censor.idx_L.begin(),censor.idx_L.end(),known),censor.idx_L.end());
  censor.idx_U.erase(std::remove_if(censor.idx_U.begin(),censor.idx_U.end(),known),censor.idx_U.end());

  // Working buffers for analytical survivor products
  censor.log_p_lower.resize(censor.n_rows);
  censor.log_p_upper.resize(censor.n_rows);

  censor.bound_work.resize(censor.n_rows);
  censor.active_rows.reserve(n_acc);
  censor.winner_row.resize(1);
  censor.breaks.reserve(n_acc+2);

  return censor;
}



