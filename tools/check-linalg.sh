#!/bin/sh
R_BIN=$1
work=$(mktemp -d "${TMPDIR:-/tmp}/emc2-linalg.XXXXXX") || exit 1
trap 'rm -rf "$work"' 0
trap 'exit 1' 1 2 15

cat > "$work/linalg.cpp" <<'EOF'
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
int main() {
  int n = 1, pivot, info;
  double x = 1;
  x = F77_CALL(ddot)(&n, &x, &n, &x, &n);
  F77_CALL(dgetrf)(&n, &n, &x, &n, &pivot, &info);
  return info;
}
EOF

link_linalg() {
  printf 'CXX_STD = CXX17\nPKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) %s\n' "$1" > "$work/Makevars"
  cat >> "$work/Makevars" <<'EOF'
all: check-link
check-link: linalg.o
	$(CXX) $(ALL_CXXFLAGS) $(LDFLAGS) -o $@ linalg.o $(ALL_LIBS)
EOF
  (cd "$work" && "$R_BIN" CMD SHLIB --preclean linalg.cpp) > "$work/link.log" 2>&1
}

if link_linalg '$(FLIBS)'; then
  echo "configure: R BLAS, LAPACK and Fortran libraries available" >&2
  printf '%s\n' '$(FLIBS)'
  exit 0
fi

# CRAN macOS R ships the runtimes needed by C/C++ packages using BLAS/LAPACK.
if [ "$(uname)" = Darwin ]; then
  r_lib="$("$R_BIN" RHOME)/lib"
  bundled=""
  for name in gfortran quadmath; do
    for runtime in "$r_lib/lib$name.dylib" "$r_lib"/lib"$name".*.dylib; do
      [ -f "$runtime" ] || continue
      bundled="$bundled \"\$(R_HOME)/lib/${runtime##*/}\""
      break
    done
  done
  if [ -n "$bundled" ] && link_linalg "$bundled"; then
    echo "configure: using R's bundled Fortran libraries" >&2
    printf '%s\n' "$bundled"
    exit 0
  fi
fi

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
