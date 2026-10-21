#ifndef EMC2_OPENMP_H
#define EMC2_OPENMP_H

#include "build_info.h"

#if defined(_OPENMP) && EMC2_OMP_PAUSE && (defined(__linux__) || defined(__APPLE__))
#include <omp.h>
#include <dlfcn.h>
#include <cstring>
#endif

namespace emc2 {

// Release GNU worker threads before returning to R, which may fork next.
struct OpenMPPoolRelease {
  ~OpenMPPoolRelease() {
#if defined(_OPENMP) && EMC2_OMP_PAUSE && (defined(__linux__) || defined(__APPLE__))
    static const bool gnu = [] {
      Dl_info info;
      if (!dladdr(reinterpret_cast<void*>(&omp_pause_resource_all), &info) || !info.dli_fname)
        return false;
      const char* name = std::strrchr(info.dli_fname, '/');
      return std::strstr(name ? name + 1 : info.dli_fname, "gomp") != nullptr;
    }();
    if (gnu) omp_pause_resource_all(omp_pause_hard);
#endif
  }
};

} // namespace emc2
#endif
