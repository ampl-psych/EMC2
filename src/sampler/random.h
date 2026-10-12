#ifndef EMC2_SAMPLER_RANDOM_H
#define EMC2_SAMPLER_RANDOM_H

#include <Rcpp.h>
#include <array>
#include <algorithm>
#include <random>
#include <cmath>
#include <sstream>

namespace emc2 {

// Seed and restore on the R thread. Each chain owns its engine and distribution state.
class Random {
  std::mt19937_64 engine;
  static std::string format() {
#if defined(_LIBCPP_VERSION)
    return "libc++ " + std::to_string(_LIBCPP_VERSION);
#elif defined(__GLIBCXX__)
    return "libstdc++ " + std::to_string(__GLIBCXX__);
#else
    return "MSVC " + std::to_string(_MSC_VER);
#endif
  }
  double spare = 0;
  bool has_spare = false;
 public:
  explicit Random(SEXP saved = R_NilValue) {
    if (!Rf_isNull(saved)) {
      Rcpp::List state(saved);
      if (state.containsElementNamed("format") && Rcpp::as<std::string>(state["format"]) == format()) {
        std::istringstream input(Rcpp::as<std::string>(state["engine"]));
        if (!(input >> engine)) Rcpp::stop("invalid sampler RNG state");
        spare = Rcpp::as<double>(state["normal"]);
        has_spare = Rcpp::as<bool>(state["has_normal"]);
        return;
      }
    }
    // A different C++ library starts a fresh stream from R's current seed.
    std::array<uint32_t, 8> words;
    for (auto& word : words) word = static_cast<uint32_t>(R::runif(0, 4294967296.0));
    std::seed_seq seed(words.begin(), words.end());
    engine.seed(seed);
  }
  Rcpp::List state() const {
    std::ostringstream output;
    output << engine;
    return Rcpp::List::create(Rcpp::Named("engine") = output.str(),
      Rcpp::Named("normal") = spare, Rcpp::Named("has_normal") = has_spare,
      Rcpp::Named("format") = format());
  }
  double rnorm(double mean, double sd) {
    if (has_spare) { has_spare = false; return mean + sd * spare; }
    // Polar Box-Muller; the spare is saved as an R double, independent of the C++ library.
    double x, y, radius;
    do {
      x = runif(-1, 1); y = runif(-1, 1);
      radius = x * x + y * y;
    } while (radius >= 1 || radius == 0);
    const double factor = std::sqrt(-2 * std::log(radius) / radius);
    spare = y * factor;
    has_spare = true;
    return mean + sd * x * factor;
  }
  double runif(double lower, double upper) {
    return lower + (upper - lower) * std::generate_canonical<double, 53>(engine);
  }
  double rgamma(double shape, double scale) {
    return std::gamma_distribution<double>(shape, scale)(engine);
  }
  double rchisq(double df) { return rgamma(df / 2, 2); }
  void rmultinom(int n, const double* p, int k, int* counts) {
    double remaining = 1;
    for (int j = 0; j + 1 < k; ++j) {
      counts[j] = n ? std::binomial_distribution<int>(n, std::min(1., p[j] / remaining))(engine) : 0;
      n -= counts[j];
      remaining -= p[j];
    }
    counts[k - 1] = n;
  }
};

} // namespace emc2
#endif
