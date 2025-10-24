# Circular Diffusion Model (CDM) - R implementation
#
# Provides:
# - rCDM: generator function to simulate RT and angle R
# - dCDM: joint density of (RT, R)
# - series_bessel_fpt: helper for the zero-drift first-passage time density
# - CDM_nll: convenience negative log-likelihood matching the notebook API
#
#
# ============================================================================

# Load built-in Bessel zeros and J1(zeros) values
#
# Returns embedded constants for zeros of J0 and J1 evaluated at those zeros.
#
load_bessel_values <- function() {
  list(zeros = c(
    2.404825557695773, 5.520078110286311, 8.653727912911013, 11.791534439014281, 14.930917708487787,
    18.071063967910924, 21.21163662987926, 24.352471530749302, 27.493479132040253, 30.634606468431976,
    33.77582021357357, 36.917098353664045, 40.05842576462824, 43.19979171317673, 46.341188371661815,
    49.482609897397815, 52.624051841115, 55.76551075501998, 58.90698392608094, 62.04846919022717,
    65.18996480020687, 68.3314693298568, 71.47298160359374, 74.61450064370183, 77.75602563038805,
    80.89755587113763, 84.0390907769382, 87.18062984364116, 90.32217263721049, 93.46371878194478,
    96.60526795099626, 99.7468198586806, 102.8883742541948, 106.02993091645162, 109.17148964980538,
    112.3130502804949, 115.45461265366694, 118.59617663087253, 121.73774208795096, 124.87930891323295,
    128.02087700600833, 131.1624462752139, 134.30401663830546, 137.44558802028428, 140.58716035285428,
    143.72873357368974, 146.87030762579664, 150.01188245695477, 153.15345801922788, 156.29503426853353,
    159.43661116426316, 162.57818866894667, 165.71976674795502, 168.86134536923583, 172.0029245030782,
    175.14450412190274, 178.28608420007376, 181.42766471373105, 184.5692456406387, 187.71082696004936,
    190.85240865258152, 193.99399070010912, 197.1355730856614, 200.2771557933324, 203.41873880819864,
    206.56032211624446, 209.70190570429406, 212.8434895599495, 215.98507367153402, 219.12665802804057,
    222.2682426190843, 225.40982743485932, 228.5514124660988, 231.69299770403853, 234.83458314038324,
    237.97616876727565, 241.11775457726802, 244.2593405632957, 247.4009267186528, 250.54251303696995,
    253.6840995121931, 256.82568613856444, 259.9672729106045, 263.1088598230955, 266.2504468710659,
    269.39203404977604, 272.5336213547049, 275.67520878153744, 278.8167963261531, 281.9583839846149,
    285.09997175315954, 288.2415596281877, 291.3831476062552, 294.524735684065, 297.66632385845895,
    300.80791212641117, 303.9495004850206, 307.09108893150506, 310.232677463195, 313.37426607752786
  ),
  JVZ = c(
    0.5191474972894667, -0.34026480655836816, 0.2714522999283819, -0.23245983136472478, 0.20654643307799606,
    -0.18772880304043937, 0.17326589422922986, -0.16170155068925005, 0.15218121377059451, -0.14416597768637318,
    0.13729694340850293, -0.13132462666866793, 0.1260694971272734, -0.12139862477175013, 0.11721119889066538,
    -0.1134291926164298, 0.10999114304627805, -0.10684788825471288, 0.10395957286936207, -0.10129349893394325,
    0.0988225538011999, -0.09652404046467991, 0.0943787939846764, -0.0923705048235533, 0.09048519416295771,
    -0.08871080244096977, 0.08703686332409762, -0.08545424291091487, 0.08395492928345759, -0.08253186130830982,
    0.08117878831953207, -0.07989015430874276, 0.07866100171930494, -0.07748689103965989, 0.07636383321829139,
    -0.07528823255205498, 0.0742568381822715, -0.07326670270620797, 0.0723151467023698, -0.071399728196232,
    0.07051821627333571, -0.06966856819003345, 0.06884890944684942, -0.068057516381685, 0.06729280091473994,
    -0.0665532971377112, 0.0658376494894271, -0.06514460230079303, 0.06447299052550669, -0.06382173150081485,
    0.06318981760571403, -0.06257630970331138, 0.06198033127024817, -0.06140106312970012, 0.060837738715961355,
    -0.06028963980834689, 0.05975609268041637, -0.05923646461756444, 0.05873016076204431, -0.058236621249651344,
    0.057755318606728306, -0.057285755379976114, 0.05682746197485655, -0.05637999468123289, 0.05594293386737833,
    -0.05551588232564261, 0.05509846375495183, -0.054690321366964154, 0.05429111660414672, -0.05390052795930567,
    0.05351824988721485, -0.053143991799969814, 0.05277747713856062, -0.05241844251392233, 0.052066636911401065,
    -0.051721820953175804, 0.05138376621371111, -0.051052254583793026, 0.05072707767912495, -0.05040803628983996,
    0.050094939867625976, -0.0497876060474634, 0.04948586020124814, -0.04918953502081812, 0.04889847012812101,
    -0.0486125117104598, 0.048331512178931856, -0.04805532984833819, 0.04778382863698613, -0.047516877784940494,
    0.04725435158939827, -0.04699612915597017, 0.04674209416475137, -0.04649213465015328, 0.046246142793549695,
    -0.04600401472786514, 0.045765650353301025, -0.045530953163457046, 0.045299830081161785, -0.04507219130337835
  ))
}

# First-passage time density for zero-drift 2D diffusion to a circular boundary
series_bessel_fpt <- function(t, a = 1, sigma = 1, bessel = NULL) {
  if (is.null(bessel)) {
    bessel <- load_bessel_values()
  }

  # Handle empty input
  if (length(t) == 0L) {
    return(numeric(0))
  }

  t <- pmax(as.numeric(t), 0)
  n <- length(t)

  # Pre-compute constants (avoiding repeated arithmetic)
  zeros <- bessel$zeros
  JVZ <- bessel$JVZ
  weights <- zeros / JVZ
  base_rates <- (zeros * zeros) / 2

  # Row-specific scaling v = sigma^2 / a^2, and s = t * v
  v <- (sigma * sigma) / (a * a)
  s <- t * v

  # exp(- s_i * base_rates_j) and then row-wise dot with weights, scaled by v_i
  exp_mat <- exp(-outer(s, base_rates, "*"))
  as.vector(v * (exp_mat %*% weights))
}

## Joint density of (RT, R) under the Circular Diffusion Model
#
# Vectorized over trials via a per-trial parameter matrix `pars` (columns:
# mu_x, mu_y, a, t0, sigma). Returns a vector of densities/log-densities of
# the same length as `rt`.
dCDM <- function(rt, R, pars) {

  n <- length(rt)
  rt    <- as.numeric(rt)
  R <- as.numeric(R)
  t0 <- pars[,"t0"]
  a <- pars[,"a"]
  mu_x <- pars[,"mu_x"]
  mu_y <- pars[,"mu_y"]
  sigma <- pars[,"sigma"]

  # First-passage time density per-trial
  bessel <- load_bessel_values()
  fpt_z <- series_bessel_fpt(rt - t0, a = a, sigma = sigma, bessel = bessel)

  cos_theta <- cos(R)
  sin_theta <- sin(R)
  mu_sq <- mu_x * mu_x + mu_y * mu_y
  term1 <- a * (mu_x * cos_theta + mu_y * sin_theta)
  term2 <- 0.5 * mu_sq * (rt - t0)

  raw_density <- exp(term1 - term2) * fpt_z
  raw_density[raw_density < 0] <- 0
  return(raw_density)
}

## Generator for the Circular Diffusion Model
#
# Simulates one (rt, R) per row of `pars` (columns: mu_x, mu_y, a, t0, sigma).
rCDM <- function(pars, ok = rep(TRUE, nrow(pars)), dt = 0.001, max_steps = 1e7L) {
  if (!is.matrix(pars)) stop("'pars' must be a matrix with columns mu_x, mu_y, a, t0, sigma")
  n <- nrow(pars)
  out <- data.frame(rt = rep(NA_real_, n), R = rep(NA_real_, n))
  if (n == 0L) return(out)
  if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0) stop("'dt' must be a single positive finite number.")

  mu_x <- as.numeric(pars[,"mu_x"]) ; mu_y <- as.numeric(pars[,"mu_y"]) ;
  a    <- as.numeric(pars[,"a"])    ; t0   <- as.numeric(pars[,"t0"])   ;
  sigma<- as.numeric(pars[,"sigma"]) ;

  idx <- which(ok)
  if (length(idx) == 0L) return(out)

  for (i in idx) {
    ax   <- a[i]; t0_i <- t0[i]; sig_i <- sigma[i]
    muxy <- c(mu_x[i], mu_y[i])
    if (!is.finite(ax) || ax <= 0) next
    if (!is.finite(t0_i) || t0_i < 0) next
    if (!is.finite(sig_i) || sig_i <= 0) next

    x <- c(0.0, 0.0)
    rt_i <- 0.0
    steps <- 0L
    sqrt_dt_sigma <- sqrt(dt) * sig_i
    drift_inc <- muxy * dt
    a_sq <- ax * ax
    repeat {
      r_sq <- x[1] * x[1] + x[2] * x[2]
      if (r_sq >= a_sq) break

      # Euler-Maruyama step: dx = mu*dt + sigma*sqrt(dt)*dW
      x <- x + drift_inc + sqrt_dt_sigma * stats::rnorm(2L)
      rt_i <- rt_i + dt
      steps <- steps + 1L
      if (steps >= max_steps) stop("Maximum number of steps reached; consider increasing 'a' or 'dt'.")
    }
    out$R[i] <- atan2(x[2], x[1])
    out$rt[i]    <- t0_i + rt_i
  }
  out
}



#' The Circular Diffusion Model (CDM)
#'
#' Model file to estimate the Circular Diffusion Model (CDM) in EMC2.
#'
#' @details
#' Parameters and defaults are on the transformed scale.
#' All parameters not specified in the `formula` of `design()` are
#' assumed constant with the defaults in `p_types`.
#'
#' | Parameter | Transform | Natural scale | Default   | Interpretation |
#' |-----------|-----------|---------------|-----------|----------------|
#' | mu_x      | identity  | (-Inf, Inf)   | 0         | Drift x-component |
#' | mu_y      | identity  | (-Inf, Inf)   | 0         | Drift y-component |
#' | a         | log       | (0, Inf)      | log(1)    | Boundary radius |
#' | t0       | log       | [0, Inf)      | log(0)    | Non-decision time |
#' | sigma     | log       | (0, Inf)      | log(1)    | Diffusion scale |
#'
#'
#' @return A model list compatible with `design()` and EMC2 fitting routines.
#' @export
CDM <- function(){
  list(
    type = "DDM",
    p_types = c(
      "mu_x" = 0,
      "mu_y" = 0,
      "a"    = log(1),
      "t0"  = log(0),
      "sigma"= log(1)
    ),
    transform = list(func = c(mu_x = "identity", mu_y = "identity",
                              a = "exp", t0 = "exp", sigma = "exp")),
    bound = list(
      minmax = cbind(
        mu_x = c(-10, 10),
        mu_y = c(-10, 10),
        a    = c(1e-4, 10),
        t0  = c(0.05, Inf),
        sigma= c(1e-4, Inf)
      )
    ),
    Ttransform = function(pars, dadm) { pars },
    # Random function: directly call rCDM with per-trial parameters
    rfun = function(data = NULL, pars) rCDM(pars, attr(pars, "ok")),
    # Density function used by the log-likelihood (captures preloaded bessel constants)
    dfun = function(rt, R, pars) dCDM(rt, R, pars),
    # DDM-style log likelihood over trials
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_ddm(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
