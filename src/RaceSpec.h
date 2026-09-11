#ifndef RACE_SPEC_H
#define RACE_SPEC_H

#include <vector>

// This is used to map parameter name to column index in a ParamTable. All default values are -1 and overwritten once created by make_race_setup (in RaceSetup)
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
  // ExG
  int col_mu = -1;
  int col_tau = -1;
  int col_sigma = -1;
};

// Choice-only models
struct ChoiceOnlySpec {
  // ParamTable column indices
  int col_location = -1;
  int col_scale    = -1;
  int col_cut      = -1;
  int col_utility  = -1;
};

// Same but for MRI models (not strictly a race model but still)
struct MRISpec {
  int col_sigma;               // last column
  int col_rho;                 // second-to-last (AR1 only)
  std::vector<int> col_means;  // all preceding columns
  int n_mean_cols;
};

// ---------------------------------------------------------------------------
// RaceScratch — flat pool of double/int slots, each of length n
// Used as a memory buffer for scatter-gather implementations
// ---------------------------------------------------------------------------
struct RaceScratch {
  std::vector<double> d;
  std::vector<int>    i;
  int n = 0;

  void reserve(int n_slots_d, int n_slots_i, int size) {
    n = size;
    d.assign(n_slots_d * size, 0.0);
    i.assign(n_slots_i * size, 0);
  }

  double* dslot(int k)       { return d.data() + k * n; }
  int*    islot(int k)       { return i.data() + k * n; }
};


// ---------------------------------------------------------------------------
// Per-model scratch slot assignments
// ---------------------------------------------------------------------------

// Shared int slots (all race models use the same layout)
namespace ScratchInt {
constexpr int idx_win   = 0;
constexpr int idx_los   = 1;
constexpr int idx_win_c = 2;   // only RDM / LBA
constexpr int idx_los_c = 3;   // only RDM / LBA
}

namespace RDM {
// double slots
constexpr int t_eff   = 0;   // noA path
constexpr int v       = 1;
constexpr int B       = 2;
constexpr int s       = 3;
constexpr int out     = 4;
constexpr int t_eff_c = 5;   // core path
constexpr int v_c     = 6;
constexpr int B_c     = 7;
constexpr int A_c     = 8;
constexpr int s_c     = 9;
constexpr int out_c   = 10;
constexpr int N_DBL   = 11;
constexpr int N_INT   = 4;   // IDX_WIN, IDX_LOS, IDX_WIN_C, IDX_LOS_C
}

namespace LBA {
// double slots — primary = core path (A > asymptotic), secondary = noA
constexpr int t_eff    = 0;
constexpr int v        = 1;
constexpr int sv       = 2;
constexpr int B        = 3;
constexpr int A        = 4;
constexpr int denom    = 5;
constexpr int out      = 6;
constexpr int t_eff_c  = 7;
constexpr int v_c      = 8;
constexpr int sv_c     = 9;   // reuses "s_c" slot semantically
constexpr int B_c      = 10;
constexpr int denom_c  = 11;
constexpr int out_c    = 12;
constexpr int N_DBL    = 13;
constexpr int N_INT    = 4;
}

namespace LNR {
// double slots — no split path
constexpr int t_eff = 0;
constexpr int m     = 1;   // lognormal mean; reuses "v" slot semantically
constexpr int s     = 2;
constexpr int out   = 3;
constexpr int N_DBL = 4;
constexpr int N_INT = 2;   // IDX_WIN, IDX_LOS only
}

namespace REXG {
// double slots — no split path
constexpr int t_eff = 0;
constexpr int mu    = 1;
constexpr int sigma = 2;
constexpr int tau   = 3;
constexpr int out   = 4;
constexpr int N_DBL = 5;
constexpr int N_INT = 2;   // IDX_WIN, IDX_LOS only
}



#endif // RACE_SPEC_H
