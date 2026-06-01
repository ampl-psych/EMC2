#ifndef RACE_SPEC_H
#define RACE_SPEC_H

#include <vector>

struct RaceSpec {
  int col_v  = -1;
  int col_B  = -1;
  int col_A  = -1;
  int col_t0 = -1;
  int col_s  = -1;
  int col_sv = -1;
  int col_m  = -1;
};

struct RaceScratch {
  std::vector<double> t_eff, v, B, A, s, sv, m, out, k, l, a;
  std::vector<int> win0, win1, los0, los1;

  void reserve(int n) {
    t_eff.resize(n); v.resize(n); B.resize(n);
    A.resize(n);     s.resize(n); out.resize(n);
    sv.resize(n);    m.resize(n); k.resize(n); l.resize(n); a.resize(n);
  }
};

#endif // RACE_SPEC_H
