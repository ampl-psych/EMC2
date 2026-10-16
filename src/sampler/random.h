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
  // The C++ library family: its distributions (gamma, binomial) are what make
  // a stream library-specific. The version is left out so that a saved object
  // continues its stream across compiler upgrades.
  static std::string format() {
#if defined(_LIBCPP_VERSION)
    return "libc++";
#elif defined(__GLIBCXX__)
    return "libstdc++";
#else
    return "MSVC";
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
    std::array<uint32_t, 8> words = r_words();
    std::seed_seq seed(words.begin(), words.end());
    engine.seed(seed);
  }
  // Eight words from R's generator; R thread only.
  static std::array<uint32_t, 8> r_words() {
    std::array<uint32_t, 8> words;
    for (auto& word : words) word = static_cast<uint32_t>(R::runif(0, 4294967296.0));
    return words;
  }
  // Reseed from the restored state and R's generator together, so that
  // set.seed() before a run picks the stream, as it did when R drew every
  // number, while a saved object still continues its own stream when the R
  // seed is the same.
  void mix(const std::array<uint32_t, 8>& words) {
    std::array<uint32_t, 16> all;
    std::copy(words.begin(), words.end(), all.begin());
    for (size_t k = 8; k < all.size(); ++k) all[k] = static_cast<uint32_t>(engine());
    std::seed_seq seed(all.begin(), all.end());
    engine.seed(seed);
    has_spare = false;
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
