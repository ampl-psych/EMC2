#ifndef RACE_SPEC_H
#define RACE_SPEC_H

#include <vector>

struct RaceSpec {
  // typically shared
  int col_v  = -1;
  int col_B  = -1;
  int col_A  = -1;
  int col_t0 = -1;
  int col_s  = -1;
  int col_sv = -1;
  // LNR
  int col_m  = -1;
  // DDM
  int col_Z  = -1;
  int col_SZ = -1;
  int col_st0= -1;
  int col_a  = -1;
};

struct RaceScratch {
  // noA gather buffers (primary — used by most entries)
  std::vector<double> t_eff, v, B, A, s, sv, m, out, k, l, a;
  // core gather buffers (secondary — used when A >= A_ASYMPTOTIC, only needed for RDM)
  std::vector<double> t_eff_c, v_c, B_c, A_c, s_c, out_c;

  std::vector<int> win0, win1, los0, los1;
  std::vector<int> idx_win0, idx_win1, idx_los0, idx_los1;
  std::vector<int> idx_win_c, idx_los_c;  // core compacted indices

  // lba precompute a denominator
  std::vector<double> denom, denom_c;

  void reserve(int n) {
    t_eff.resize(n); v.resize(n); B.resize(n);
    A.resize(n);     s.resize(n); out.resize(n);
    sv.resize(n);    m.resize(n); k.resize(n); l.resize(n); a.resize(n);
    t_eff_c.resize(n); v_c.resize(n); B_c.resize(n);
    A_c.resize(n);     s_c.resize(n); out_c.resize(n);
    win0.resize(n);  win1.resize(n);
    los0.resize(n);  los1.resize(n);
    idx_win0.resize(n); idx_win1.resize(n);
    idx_los0.resize(n); idx_los1.resize(n);
    idx_win_c.resize(n); idx_los_c.resize(n);
    denom.resize(n); denom_c.resize(n);
  }
};

#endif // RACE_SPEC_H
