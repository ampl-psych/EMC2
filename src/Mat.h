#ifndef MAT_H
#define MAT_H

#include <Rcpp.h>          // boundary only — not used in hot paths


// =============================================================================
// Mat  —  plain column-major matrix; no Rcpp dependency
// =============================================================================

struct Mat {
  int nrow = 0, ncol = 0;
  std::vector<double> data;   // column-major: element (r,c) = data[c*nrow + r]

  Mat() = default;
  Mat(int r, int c) : nrow(r), ncol(c), data(r * c, 0.0) {}

  double*       colptr(int j)       { return data.data() + j * nrow; }
  const double* colptr(int j) const { return data.data() + j * nrow; }

  double&       operator()(int r, int c)       { return data[c * nrow + r]; }
  const double& operator()(int r, int c) const { return data[c * nrow + r]; }

  // Construct from a bare SEXP (real matrix) — the only Rcpp-free R boundary
  static Mat from_sexp(SEXP s) {
    if (!Rf_isReal(s))
      Rf_error("Mat::from_sexp: expected a numeric matrix");
    int nr = Rf_nrows(s), nc = Rf_ncols(s);
    Mat m(nr, nc);
    const double* src = REAL(s);
    std::copy(src, src + nr * nc, m.data.data());
    return m;
  }
};

#endif // mat_h
