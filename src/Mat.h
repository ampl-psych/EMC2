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

  // Construct from Rcpp::NumericMatrix (at the R boundary, before parallel region)
  static Mat from_rcpp(const Rcpp::NumericMatrix& m) {
    return from_sexp(m);
  }

  Mat clone() const {
    return *this;  // std::vector<double> data is deep-copied by value
  }
};


// Boolean version, use uint8_t to allow for raw pointer access (bool doesnt do that for some reason)
struct MatBool {
  int nrow = 0, ncol = 0;
  std::vector<uint8_t> data;        // owning storage — empty when view
  const uint8_t* view_ptr = nullptr; // non-owning — set when constructed as view

  MatBool() = default;
  MatBool(int r, int c, uint8_t fill = 1)
    : nrow(r), ncol(c), data(r * c, fill), view_ptr(nullptr) {}

  // Non-owning view into a single column of another MatBool.
  // Lifetime: caller must ensure the source outlives this view.
  static MatBool col_view(const MatBool& src, int col) {
    MatBool v;
    v.nrow     = src.nrow;
    v.ncol     = 1;
    v.view_ptr = src.colptr(col);
    return v;
  }

  uint8_t*       colptr(int j)       { return (view_ptr ? const_cast<uint8_t*>(view_ptr) : data.data()) + j * nrow; }
  const uint8_t* colptr(int j) const { return (view_ptr ? view_ptr : data.data()) + j * nrow; }

  uint8_t&       operator()(int r, int c)       { return colptr(c)[r]; }
  const uint8_t& operator()(int r, int c) const { return colptr(c)[r]; }

  MatBool clone() const {
    if (view_ptr) {
      // materialise the view into an owned copy
      MatBool m(nrow, ncol);
      std::copy(view_ptr, view_ptr + nrow * ncol, m.data.data());
      return m;
    }
    return *this;
  }
};



#endif // mat_h
