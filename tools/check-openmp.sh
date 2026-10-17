# Sourced by configure; use the same C++17 compiler and user flags as the package.
OMP_CXXFLAGS=""
OMP_LIBS=""
OMP_PAUSE=0
OMP_TMP=$(mktemp -d "${TMPDIR:-/tmp}/emc2-openmp.XXXXXX") || exit 1
trap 'rm -rf "$OMP_TMP"' 0
trap 'exit 1' 1 2 15

cat > "$OMP_TMP/conftest.cpp" <<'EOF'
#include <omp.h>
#if !defined(_OPENMP) || _OPENMP < 201307
#error OpenMP 4.0 is required for task groups
#endif
int main() {
  int value = 0;
  #pragma omp parallel
  #pragma omp single
  #pragma omp taskgroup
  {
    #pragma omp task shared(value)
    { value = 1; }
  }
  return value != 1 || omp_get_max_threads() < 1;
}
EOF

check_openmp() {
  printf 'CXX_STD = CXX17\nPKG_CXXFLAGS = %s %s\nPKG_LIBS = %s\n' \
    "$1" "$EXTRA_CXXFLAGS" "$2" > "$OMP_TMP/Makevars"
  cat >> "$OMP_TMP/Makevars" <<'EOF'
all: check-link
check-link: conftest.o
	$(CXX) $(ALL_CXXFLAGS) $(LDFLAGS) -o $@ conftest.o $(ALL_LIBS)
EOF
  (cd "$OMP_TMP" && "$R_BIN" CMD SHLIB --preclean conftest.cpp) >> config.log 2>&1
}

if [ "$ENABLE_OPENMP" = yes ]; then
  if check_openmp '$(SHLIB_OPENMP_CXXFLAGS)' '$(SHLIB_OPENMP_CXXFLAGS)'; then
    OMP_CXXFLAGS='$(SHLIB_OPENMP_CXXFLAGS)'
    OMP_LIBS='$(SHLIB_OPENMP_CXXFLAGS)'
  elif [ "$(uname)" = Darwin ]; then
    # Prefer R's runtime; standard R build tools and Homebrew supply omp.h.
    for prefix in /usr/local /opt/homebrew/opt/libomp /usr/local/opt/libomp; do
      flags="-Xclang -fopenmp -I${prefix}/include"
      libs="-L\"\$(R_HOME)/lib\" -L${prefix}/lib -lomp"
      if check_openmp "$flags" "$libs"; then
        OMP_CXXFLAGS=$flags
        OMP_LIBS=$libs
        break
      fi
    done
  fi
fi

if [ -n "$OMP_LIBS" ]; then
  # GCC 9 provides the pool-release API but reports OpenMP 4.5.
  cat > "$OMP_TMP/conftest.cpp" <<'EOF'
#include <omp.h>
int main() { return omp_pause_resource_all(omp_pause_hard); }
EOF
  if check_openmp "$OMP_CXXFLAGS" "$OMP_LIBS"; then OMP_PAUSE=1; fi
  echo "configure: OpenMP enabled"
else
  echo "configure: building without OpenMP (serial sampler)"
fi
