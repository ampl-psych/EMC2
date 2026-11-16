# Independent formulas using R's normal probabilities, not EMC2's ex-Gaussian
# helpers. Separate intervals resolve the Gaussian rise and exponential tail.
exg_integral_reference <- function(p, survivor, upper = Inf) {
  mu <- p["muS"]; sigma <- p["sigmaS"]; tau <- p["tauS"]
  lower <- p["exgS_lb"]
  if (upper <= lower) return(0)
  raw_survivor <- function(t) {
    if (t == -Inf) return(1)
    z <- (t - mu) / sigma
    pnorm(-z) + exp(sigma / tau * (sigma / (2 * tau) - z) +
                       pnorm(z - sigma / tau, log.p = TRUE))
  }
  normalizer <- raw_survivor(lower)
  integrand <- function(t) {
    z <- (t - mu) / sigma
    exp(sigma / tau * (sigma / (2 * tau) - z) +
          pnorm(z - sigma / tau, log.p = TRUE)) / tau / normalizer *
      survivor(t + p["SSD"])
  }
  cuts <- c(mu + sigma * c(-12, -4, 0, 4, 12), mu + tau * c(1, 5, 20))
  cuts <- sort(unique(c(lower, cuts[cuts > lower & cuts < upper], upper)))
  sum(vapply(seq_len(length(cuts) - 1L), function(i)
    integrate(integrand, cuts[i], cuts[i + 1L], abs.tol = 1e-12,
              rel.tol = 1e-10, subdivisions = 500L)$value, 0.0))
}

test_that("ex-Gaussian stop integrals resolve the saved-draw narrow transitions", {
  # Participant 1, draw 1500: chains 1/3, trials 2/24. No saved fit is needed.
  cases <- list(
    list(mu = c(.496254136774276, .555244994937601),
         sigma = c(.0661910153761801, .11809735733019),
         tau = c(.474829368837576, 2.26696732580507),
         muS = .196294421008146, sigmaS = .000139180792335576,
         tauS = .052014549083965, SSD = .2),
    list(mu = c(.502078583481696, .717034489190613),
         sigma = c(.0658184098505406, .201023385230101),
         tau = c(.431362489568333, 1.61749665466957),
         muS = .187178522270336, sigmaS = .000125385436673718,
         tauS = .0464182481895738, SSD = .05)
  )
  for (case in cases) {
    pars <- as.matrix(as.data.frame(case))
    pars <- cbind(pars, exg_lb = 0, exgS_lb = 0,
                  exg_normal_ratio = 0, exg_shifted_exp_sigma = 0)
    survivor <- function(t) {
      result <- rep(1, length(t))
      for (i in 1:2) {
        mu <- case$mu[i]; sig <- case$sigma[i]; tau <- case$tau[i]
        raw <- function(t) {
          z <- (t - mu) / sig
          pnorm(-z) + exp(sig / tau * (sig / (2 * tau) - z) +
                             pnorm(z - sig / tau, log.p = TRUE))
        }
        result <- result * raw(pmax(t, 0)) / raw(0)
      }
      result
    }
    # Includes both sides of the former likelihood discontinuity and an
    # ordinary sigma. The fixed window omits at most exp(-16) stop mass here.
    for (sigma in c(case$sigmaS, .0001, .00014, .00015, .001, .03)) {
      pars[, "sigmaS"] <- sigma
      expected <- exg_integral_reference(pars[1, ], survivor)
      actual <- ss_stop_probability(pars, 2L, "SSEXG", Inf)
      expect_true(actual >= 0 && actual <= 1)
      expect_lt(abs(actual - expected), 1.2e-7)
    }
  }
})

test_that("narrow ex-Gaussian stop integrals retain truncation and deadlines", {
  for (go in c("SSEXG", "SSRDEX")) {
    if (go == "SSEXG") {
      g <- c(mu = .5, sigma = .1, tau = .1, exg_lb = 0)
      # Independent ex-Gaussian survivor for these ordinary go parameters.
      raw <- function(t) pnorm((.5 - t) / .1) +
        exp(.5 - (t - .5) / .1 + pnorm((t - .5) / .1 - 1, log.p = TRUE))
      survivor <- function(t) (raw(pmax(t, 0)) / raw(0))^2
    } else {
      g <- c(v = 3, B = 1.2, A = 0, t0 = .15, s = 1)
      survivor <- function(t) {
        dt <- pmax(t - .15, .Machine$double.eps)
        (pnorm((1.2 - 3 * dt) / sqrt(dt)) -
           exp(2 * 3 * 1.2) * pnorm(-(1.2 + 3 * dt) / sqrt(dt)))^2
      }
    }
    for (lower in c(-Inf, 0, .2, .21)) {
      p <- c(g, muS = .2, sigmaS = .00013, tauS = .05, exgS_lb = lower,
             exg_normal_ratio = 0, exg_shifted_exp_sigma = 0, SSD = .2)
      pars <- rbind(p, p)
      # Before, inside and after the density transition; also beyond it when
      # the truncation lower bound removes the transition entirely.
      for (upper in c(0, .1999, .2, .2001, .22, Inf)) {
        actual <- ss_stop_probability(pars, 2L, go, upper)
        expected <- exg_integral_reference(p, survivor, upper)
        expect_true(abs(actual - expected) < 1.2e-7,
                    info = paste(go, lower, upper))
      }
    }
  }
})
