#!/bin/sh
R_BIN=$1
work=$(mktemp -d "${TMPDIR:-/tmp}/emc2-linalg.XXXXXX") || exit 1
trap 'rm -rf "$work"' 0

cat > "$work/Makevars" <<'EOF'
CXX_STD = CXX17
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
EOF
cat > "$work/linalg.cpp" <<'EOF'
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
extern "C" void emc2_linalg_check() {
  int n = 1, pivot, info;
  double x = 1;
  x = F77_CALL(ddot)(&n, &x, &n, &x, &n);
  F77_CALL(dgetrf)(&n, &n, &x, &n, &pivot, &info);
}
EOF

if (cd "$work" && "$R_BIN" CMD SHLIB linalg.cpp) > "$work/link.log" 2>&1; then
  echo "configure: R BLAS, LAPACK and Fortran libraries available"
else
  cat "$work/link.log" >&2
  cat >&2 <<'EOF'
configure: error: cannot link against R's BLAS, LAPACK and Fortran libraries.
Install the build tools matching this R installation and retry:
  macOS: Xcode Command Line Tools and Fortran from https://mac.r-project.org/tools/
  Windows: matching Rtools from https://cran.r-project.org/bin/windows/Rtools/
  Linux/Unix: R development packages and the required compiler runtime libraries.
See the compiler/linker output above for the missing or incompatible dependency.
EOF
  exit 1
fi
