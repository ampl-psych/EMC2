// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//--------------------------------------------------------------------
// Minimal Hungarian assignment solver for small square cost matrices
//--------------------------------------------------------------------
static uvec hungarian(const mat& cost)
{
  const unsigned int n = cost.n_rows;              // (square!)
  const double INF = std::numeric_limits<double>::infinity();

  vec u(n + 1, fill::zeros), v(n + 1, fill::zeros);
  uvec p(n + 1, fill::zeros), way(n + 1, fill::zeros);

  for (unsigned int i = 1; i <= n; ++i) {
    p(0) = i;
    unsigned int j0 = 0;                            // current column
    vec minv(n + 1, fill::value(INF));
    uvec used(n + 1, fill::zeros);

    do {
      used(j0) = 1;
      unsigned int i0 = p(j0),  j1 = 0;
      double delta = INF;

      for (unsigned int j = 1; j <= n; ++j) if (!used(j)) {
        double cur = cost(i0 - 1, j - 1) - u(i0) - v(j);
        if (cur < minv(j)) { minv(j) = cur; way(j) = j0; }
        if (minv(j) < delta) { delta = minv(j); j1 = j; }
      }

      for (unsigned int j = 0; j <= n; ++j) {
        if (used(j)) { u(p(j)) += delta; v(j) -= delta; }
        else         { minv(j)  -= delta; }
      }
      j0 = j1;
    } while (p(j0) != 0);

    // augmenting path
    do {
      unsigned int j1 = way(j0);
      p(j0) = p(j1);
      j0 = j1;
    } while (j0 != 0);
  }

  uvec assign(n);
  for (unsigned int j = 1; j <= n; ++j) {
    unsigned int i = p(j);
    assign(i - 1) = j - 1;   // 0‑based
  }
  return assign;
}

//--------------------------------------------------------------------
// sp_new_cpp : drop‑in replacement for original R sp_new() helper
//--------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List sp_new(const int        iter,            // 1‑based MCMC slice index
                      const arma::cube& lambda_varimax, // p × q × M
                      const int        q,               // # factors (cols)
                      const int        p,               // # variables (rows)
                      const int        dim_all_c,       // 2^q   (rows of all_c)
                      const arma::mat& all_c,           // dim_all_c × q  sign grid
                      const arma::mat& lambda_hat,      // p × q centroid
                      const arma::uvec& st,             // *unused* (kept for sig.)
                      arma::mat       cost_matrix,      // q × q scratch  (ignored)
                      arma::mat       perm)             // dim_all_c × q scratch (ignored)
{
  // -----------------------------------------------------------------
  // Extract current loading matrix
  mat lambda = lambda_varimax.slice(iter - 1);          // 0‑based slice

  double best_cost = std::numeric_limits<double>::infinity();
  rowvec best_c(q, fill::ones);
  uvec   best_perm(q);
  mat    best_switch(p, q, fill::zeros);

  // work buffers
  mat switch_centroid(p, q);
  mat local_cost(q, q);

  // -----------------------------------------------------------------
  for (int idx = 0; idx < dim_all_c; ++idx) {
    rowvec c_vec = all_c.row(idx);

    // apply sign flip to centroid: column‑wise multiplication
    switch_centroid = lambda_hat;                      // copy
    for (int j = 0; j < q; ++j) switch_centroid.col(j) *= c_vec(j);

    // build cost matrix
    for (int j = 0; j < q; ++j) {
      for (int k = 0; k < q; ++k) {
        local_cost(j, k) = accu(square(lambda.col(k) - switch_centroid.col(j)));
      }
    }

    // assignment by Hungarian algorithm
    uvec perm_vec = hungarian(local_cost);
    double cost = 0.0;
    for (int j = 0; j < q; ++j) cost += local_cost(j, perm_vec(j));

    // retain best solution
    if (cost < best_cost) {
      best_cost  = cost;
      best_c     = c_vec;
      best_perm  = perm_vec;
    }
  }

  // Build switchedMatrix for best solution
  for (int j = 0; j < q; ++j) {
    best_switch.col(j) = best_c(j) * lambda.col(best_perm(j));
  }

  // Convert permutation to 1‑based for R side
  IntegerVector min_perm(q);
  for (int j = 0; j < q; ++j) min_perm[j] = best_perm(j) + 1;

  return List::create(
    _["min_perm"]       = min_perm,
    _["min_c"]          = NumericVector(best_c.begin(), best_c.end()),
    _["min_cost"]       = best_cost,
    _["switchedMatrix"] = best_switch
  );
}
