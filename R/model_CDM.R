# Circular Diffusion Model (CDM)
#
# Load built-in Bessel zeros and J1(zeros) values
# Returns embedded constants for zeros of J0 and J1 evaluated at those zeros.
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
# --- First-passage time density for zero-drift 2D diffusion to a circular boundary
series_bessel_fpt <- function(t, a = 1, sigma = 1, bessel = NULL) {
  if (is.null(bessel)) bessel <- load_bessel_values()
  if (length(t) == 0L) return(numeric(0))

  t <- pmax(as.numeric(t), 0)

  zeros   <- bessel$zeros          # zeros of J0
  JVZ     <- bessel$JVZ            # J1(zeros)
  lambda  <- (zeros * zeros) / 2   # eigenvalues in s = t * sigma^2 / a^2 units

  # Raw spectral weights from J0/J1 structure
  weights <- zeros / JVZ
  # norm_const <- sum(raw_weights / lambda)
  # weights    <- raw_weights / norm_const

  # Per-trial scaling: s = t * (sigma^2 / a^2)
  v_scale <- (sigma * sigma) / (a * a)
  s <- t * v_scale

  # Classical eigen-sum for FPT density (zero drift, marginal over angle):
  # f0(t) = v_scale * sum_n w_n exp(-lambda_n s)
  exp_mat <- exp(-outer(s, lambda, "*"))
  out <- as.vector(v_scale * (exp_mat %*% weights))

  # Numerical safety: clip tiny negative values from truncation/roundoff
  out[out < 0] <- 0
  out
}

# --- Small-t safeguard approximation for zero-drift FPT (unchanged, used only in a very small window)
small_t_fpt <- function(t_scaled, x_scaled, bessel = NULL) {
  if (is.null(bessel)) {
    bessel <- load_bessel_values()
  }
  zeros1 <- bessel$zeros[1]
  t <- pmax(as.numeric(t_scaled), 0)
  x <- as.numeric(x_scaled)

  term1 <- ((1 - x) * (1 + t)^2) / (sqrt(x + t) * (t^(1.5)))
  term2 <- exp(-0.5 * (1 - x)^2 / t - 0.5 * zeros1 * zeros1 * t)
  out <- term1 * term2
  out[!is.finite(out)] <- 0
  out
}

## Joint density of (RT, R) under the Circular Diffusion Model
dCDM <- function(rt, R, pars) {

  n  <- length(rt)
  rt <- as.numeric(rt)
  R  <- as.numeric(R)

  t0    <- pars[, "t0"]
  a     <- pars[, "a"]
  v     <- pars[, "v"]
  theta <- pars[, "theta"]
  # map theta from (0,1) to [-pi, pi]
  theta <- (theta - 0.5) * 2 * pi

  sigma <- pars[, "s"]
  sv    <- if ("sv" %in% colnames(pars)) as.numeric(pars[, "sv"]) else rep(0, n)

  # Zero-drift FPT with robust small-s switch
  bessel <- load_bessel_values()
  tt <- pmax(rt - t0, 0)

  # scale s = t * sigma^2 / a^2; smooth small-t handling using blending window
  v_scale <- (sigma * sigma) / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02

  fpt_series <- series_bessel_fpt(tt, a = a, sigma = sigma, bessel = bessel)
  fpt_small  <- v_scale * small_t_fpt(
    t_scaled = s,
    x_scaled = 1e-8 / (a * a),
    bessel = bessel
  )

  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  fpt <- (1 - w) * fpt_small + w * fpt_series

  # safety: nonnegative and finite
  fpt[!is.finite(fpt) | (fpt < 0)] <- 0
  log_fpt <- log(pmax(fpt, .Machine$double.xmin))

  # Angle terms
  cosR <- cos(R)
  sinR <- sin(R)
  cosT <- cos(theta)
  sinT <- sin(theta)
  cosRth <- cosR * cosT + sinR * sinT

  v2 <- v * v
  log_raw <- rep(-Inf, n)

  # sv == 0 : constant-drift likelihood (Girsanov at stopping time)
  idx0 <- (sv <= 0)
  if (any(idx0)) {
    sig2_0    <- sigma[idx0]^2
    tt0       <- tt[idx0]
    log_time0 <- -0.5 * (v2[idx0] * tt0) / sig2_0
    ang0      <- (a[idx0] * v[idx0] * cosRth[idx0]) / sig2_0
    log_raw[idx0] <- log_fpt[idx0] + log_time0 + ang0 - log(2 * pi)
  }

  # sv > 0 : integrated Gaussian drift variability
  idxv <- (sv > 0)
  if (any(idxv)) {
    tt_v <- pmax(tt[idxv], .Machine$double.eps)
    sig2 <- sigma[idxv]^2
    sv2  <- sv[idxv]^2
    D    <- sig2 + sv2 * tt_v

    # Determinant term for integrating over Gaussian drift variability.
    # The prior contributes 1/(2*pi*sv^2); after combining with the quadratic
    # form we get a factor proportional to sig2 / D (not sig2 * sv2 / D).
    log_det <- -log(D) + log(sig2)
    # radius, angular, time terms
    log_rad <- 0.5 * (a[idxv]^2 * sv2) / (sig2 * D)
    log_ang <- (a[idxv] * v[idxv] * cosRth[idxv]) / D
    log_tim <- -0.5 * (v2[idxv] * tt_v) / D

    log_raw[idxv] <- log_fpt[idxv] + log_det + log_rad + log_ang + log_tim - log(2 * pi)
  }

  dens <- exp(log_raw)
  dens[!is.finite(dens)] <- 0
  dens
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
#' | v         | log       | [0, Inf)      | log(1)    | Drift magnitude |
#' | theta     | pnorm     | [0, 1)        | qnorm(.5) | Drift direction, mapped to -pi, pi |
#' | a         | log       | [0, Inf)      | log(1)    | Boundary radius |
#' | t0        | log       | [0, Inf)      | log(0)    | Non-decision time |
#' | s         | log       | [0, Inf)      | log(1)    | Diffusion scale |
#' | sv        | log       | [0, Inf)      | log(0)    | Drift SD across trials |
#'
#'
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @export
CDM <- function(){
  list(
    type = "CDM",
    c_name = "CDM",
    p_types = c(
      "v"     = log(1),
      "theta" = qnorm(0.5),
      "a"    = log(1),
      "t0"  = log(0),
      "s"= log(1),
      "sv"   = log(0)
    ),
    transform = list(func = c(v = "exp", theta = "pnorm",
                              a = "exp", t0 = "exp", s = "exp", sv = "exp")),
    bound = list(
      minmax = cbind(
        v     = c(0, 10),
        theta = c(1e-4, 1-1e-4),
        a    = c(1e-4, 7),
        t0  = c(0.05, Inf),
        s= c(1e-4, Inf),
        sv   = c(0, 2)
      ),
      exception = c(sv = 0, t0 = 0)
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

# --- Embedded constants: zeros of J1 and J2(zeros) for HSDM FPT expansion
load_bessel_values_hsdm <- function() {
  list(
    zeros = c(
      3.831705970207512, 7.015586669815619, 10.173468135062722, 13.323691936314223, 16.470630050877634,
      19.615858510468243, 22.760084380592772, 25.903672087618382, 29.046828534916855, 32.189679910974405,
      35.332307550083868, 38.474766234771614, 41.617094212814450, 44.759318997652819, 47.901460887185436,
      51.043535183571500, 54.185553641061318, 57.327525437901009, 60.469457845347492, 63.611356698481231,
      66.753226734098490, 69.895071837495777, 73.036895225573829, 76.178699584641464, 79.320487175476302,
      82.462259914373561, 85.604019436350228, 88.745767144926305, 91.887504251694992, 95.029231808044699,
      98.170950730790778, 101.312661823038724, 104.454365791282754, 107.596063259509165, 110.737754780899209,
      113.879440847594992, 117.021121898892432, 120.162798328149009, 123.304470488635715, 126.446138698516557,
      129.587803245104025, 132.729464388509027, 135.871122364789706, 139.012777388659985, 142.154429656213838,
      145.296079347706266, 148.437726629835867, 151.579371657401489, 154.721014574856434, 157.862655517336740,
      161.004294611826280, 164.145931978044970, 167.287567729224087, 170.429201972793956, 173.570834810998010,
      176.712466341441674, 179.854096657583579, 182.995725849174829, 186.137354002649290, 189.278981201468048,
      192.420607526420260, 195.562233055884391, 198.703857865059140, 201.845482026168950, 204.987105608647740,
      208.128728679302836, 211.270351302463352, 214.411973540114992, 217.553595451992613, 220.695217095640451,
      223.836838526444509, 226.978459797641182, 230.120080960307605, 233.261702063332424, 236.403323153369491,
      239.544944274775897, 242.686565469534848, 245.828186777170782, 248.969808234658737, 252.111429876328191,
      255.253051733763616, 258.394673835707063, 261.536296208956958, 264.677918878251850, 267.819541866141725,
      270.961165193390761, 274.102788879980951, 277.244412943796139, 280.386037400863942, 283.527662265030063,
      286.669287548764783, 289.810913263240058, 292.952539418543975, 296.094166023740835, 299.235793087051500,
      302.377420615913833, 305.519048617034536, 308.660677096436061, 311.802306059500997, 314.943935510014450
    ),
    JVZ = c(
      0.402759395702553, -0.300115752526133, 0.249704877057843, -0.218359407247873, 0.196465371468657,
      -0.180063375344316, 0.167184600473818, -0.156724986252852, 0.148011109972778, -0.140605798183982,
      0.134211240310001, -0.128616622072070, 0.123667960769837, -0.119249812010690, 0.115273694120168,
      -0.111670496859211, 0.108385348943683, -0.105374055395235, 0.102600567103397, -0.100035146811523,
      0.097653015783173, -0.095433339020535, 0.093358453290455, -0.091413272155921, 0.089584821964856,
      -0.087861876039410, 0.086234663413288, -0.084694634803424, 0.083234272982226, -0.081846937926486,
      0.080526739448403, -0.079268431724519, 0.078067325407949, -0.076919213961391, 0.075820311569167,
      -0.074767200537075, 0.073756786512857, -0.072786260189239, 0.071853064408847, -0.070954865793037,
      0.070089530177261, -0.069255101263766, 0.068449782005188, -0.067671918315546, 0.066919984772397,
      -0.066192572028753, 0.065488375698257, -0.064806186514098, 0.064144881622988, -0.063503416886855,
      0.062880819010306, -0.062276178395135, 0.061688642609329, -0.061117410399967, 0.060561726189331,
      -0.060020875003365, 0.059494178619212, -0.058980992688552, 0.058480704121110, -0.057992728618295,
      0.057516508333579, -0.057051509641129, 0.056597221000638, -0.056153151900904, 0.055718831878490,
      -0.055293809601609, 0.054877651997685, -0.054469943417334, 0.054070284830448, -0.053678292046322,
      0.053293594013101, -0.052915831191146, 0.052544654000404, -0.052179721334510, 0.051820699139213,
      -0.051467259047462, 0.051119077066666, -0.050775832320041, 0.050437205836130, -0.050102879375835,
      0.049772534296667, -0.049445850453304, 0.049122504132751, -0.048802166022653, 0.048484499208068,
      -0.048169157195273, 0.047855782959665, -0.047543007989501, 0.047230451324987, -0.046917718593741,
      0.046604400039926, -0.046290068552468, 0.045974277690954, -0.045656560605942, 0.045336429769369,
      -0.045013375412588, 0.044686864603526, -0.044356340032976, 0.044021218350483, -0.043680888015461
    )
  )
}

# --- Fixed-boundary zero-drift FPT for 3D spherical diffusion (SDM)
series_sdm_fpt <- function(t, a = 1, sigma = 1, max_n = 500L) {
  if (length(t) == 0L) return(numeric(0))
  t <- pmax(as.numeric(t), 0)
  n <- seq_len(max_n)
  lambda <- 0.5 * pi * pi * (n * n)
  weights <- ((-1)^(n + 1L)) * pi * pi * (n * n)
  v_scale <- (sigma * sigma) / (a * a)
  s <- t * v_scale
  out <- as.vector(v_scale * (exp(-outer(s, lambda, "*")) %*% weights))
  out[out < 0] <- 0
  out
}

small_t_fpt_sdm <- function(t_scaled, x_scaled) {
  t <- pmax(as.numeric(t_scaled), 0)
  x <- as.numeric(x_scaled)
  term1 <- ((1 - x) * (1 + t)^2.5) / ((x + t) * (t^1.5))
  term2 <- exp(-0.5 * (1 - x)^2 / t - 0.5 * pi^2 * t)
  out <- term1 * term2
  out[!is.finite(out)] <- 0
  out
}

# --- Fixed-boundary zero-drift FPT for 4D hyperspherical diffusion (HSDM)
series_hsdm_fpt <- function(t, a = 1, sigma = 1, bessel = NULL) {
  if (is.null(bessel)) bessel <- load_bessel_values_hsdm()
  if (length(t) == 0L) return(numeric(0))
  t <- pmax(as.numeric(t), 0)
  zeros <- bessel$zeros
  JVZ <- bessel$JVZ
  lambda <- (zeros * zeros) / 2
  weights <- (zeros * zeros) / JVZ
  v_scale <- (sigma * sigma) / (a * a)
  s <- t * v_scale
  out <- as.vector(v_scale * (exp(-outer(s, lambda, "*")) %*% weights))
  out[out < 0] <- 0
  out
}

small_t_fpt_hsdm <- function(t_scaled, x_scaled, bessel = NULL) {
  if (is.null(bessel)) bessel <- load_bessel_values_hsdm()
  z1 <- bessel$zeros[1]
  t <- pmax(as.numeric(t_scaled), 0)
  x <- as.numeric(x_scaled)
  term1 <- ((1 - x) * (1 + t)^3) / ((x + t) * sqrt(x + t) * (t^1.5))
  term2 <- exp(-0.5 * (1 - x)^2 / t - 0.5 * z1 * z1 * t)
  out <- term1 * term2
  out[!is.finite(out)] <- 0
  out
}

cdm_log_density_core <- function(log_fpt, tt, a, v, sigma, sv, cos_term, dim, log_surface) {
  n <- length(tt)
  log_raw <- rep(-Inf, n)
  sig2 <- sigma * sigma
  v2 <- v * v

  idx0 <- (sv <= 0)
  if (any(idx0)) {
    log_raw[idx0] <- log_fpt[idx0] -
      0.5 * (v2[idx0] * tt[idx0]) / sig2[idx0] +
      (a[idx0] * v[idx0] * cos_term[idx0]) / sig2[idx0] -
      log_surface
  }

  idxv <- (sv > 0)
  if (any(idxv)) {
    tt_v <- pmax(tt[idxv], .Machine$double.eps)
    sig2_v <- sig2[idxv]
    sv2 <- sv[idxv]^2
    D <- sig2_v + sv2 * tt_v
    log_det <- 0.5 * dim * (log(sig2_v) - log(D))
    log_rad <- 0.5 * (a[idxv]^2 * sv2) / (sig2_v * D)
    log_ang <- (a[idxv] * v[idxv] * cos_term[idxv]) / D
    log_tim <- -0.5 * (v2[idxv] * tt_v) / D
    log_raw[idxv] <- log_fpt[idxv] + log_det + log_rad + log_ang + log_tim - log_surface
  }

  log_raw
}

log_i0_stable <- function(x) {
  ax <- abs(x)
  ax + log(pmax(besselI(ax, nu = 0, expon.scaled = TRUE), .Machine$double.xmin))
}

# Joint density of (RT, R, R2) under the Spherical Diffusion Model
dSDM <- function(rt, R, R2, pars) {
  n <- length(rt)
  rt <- as.numeric(rt)
  R <- as.numeric(R)
  R2 <- as.numeric(R2)

  t0 <- pars[, "t0"]
  a <- pars[, "a"]
  v <- pars[, "v"]
  theta1 <- pars[, "theta1"] * pi
  theta2 <- (pars[, "theta2"] - 0.5) * 2 * pi
  sigma <- pars[, "s"]
  sv <- if ("sv" %in% colnames(pars)) as.numeric(pars[, "sv"]) else rep(0, n)

  tt <- pmax(rt - t0, 0)
  v_scale <- (sigma * sigma) / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02

  fpt_series <- series_sdm_fpt(tt, a = a, sigma = sigma)
  fpt_small <- v_scale * small_t_fpt_sdm(s, 1e-8 / (a * a))
  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  fpt <- (1 - w) * fpt_small + w * fpt_series
  fpt[!is.finite(fpt) | (fpt < 0)] <- 0
  log_fpt <- log(pmax(fpt, .Machine$double.xmin))

  x0 <- cos(R)
  x1 <- sin(R) * cos(R2)
  x2 <- sin(R) * sin(R2)

  m0 <- cos(theta1)
  m1 <- sin(theta1) * cos(theta2)
  m2 <- sin(theta1) * sin(theta2)

  cos_term <- x0 * m0 + x1 * m1 + x2 * m2
  log_raw <- cdm_log_density_core(
    log_fpt = log_fpt,
    tt = tt,
    a = a,
    v = v,
    sigma = sigma,
    sv = sv,
    cos_term = cos_term,
    dim = 3,
    log_surface = log(4 * pi)
  )

  dens <- exp(log_raw)
  dens[!is.finite(dens)] <- 0
  dens
}

# Joint density of (RT, R, R2, R3) under the Hyper-Spherical Diffusion Model
dHSDM <- function(rt, R, R2, R3, pars) {
  n <- length(rt)
  rt <- as.numeric(rt)
  R <- as.numeric(R)
  R2 <- as.numeric(R2)
  R3 <- as.numeric(R3)

  t0 <- pars[, "t0"]
  a <- pars[, "a"]
  v <- pars[, "v"]
  theta1 <- pars[, "theta1"] * pi
  theta2 <- pars[, "theta2"] * pi
  theta3 <- (pars[, "theta3"] - 0.5) * 2 * pi
  sigma <- pars[, "s"]
  sv <- if ("sv" %in% colnames(pars)) as.numeric(pars[, "sv"]) else rep(0, n)

  tt <- pmax(rt - t0, 0)
  v_scale <- (sigma * sigma) / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02

  bessel <- load_bessel_values_hsdm()
  fpt_series <- series_hsdm_fpt(tt, a = a, sigma = sigma, bessel = bessel)
  fpt_small <- v_scale * small_t_fpt_hsdm(s, 1e-8 / (a * a), bessel = bessel)
  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  fpt <- (1 - w) * fpt_small + w * fpt_series
  fpt[!is.finite(fpt) | (fpt < 0)] <- 0
  log_fpt <- log(pmax(fpt, .Machine$double.xmin))

  x0 <- cos(R)
  x1 <- sin(R) * cos(R2)
  x2 <- sin(R) * sin(R2) * cos(R3)
  x3 <- sin(R) * sin(R2) * sin(R3)

  m0 <- cos(theta1)
  m1 <- sin(theta1) * cos(theta2)
  m2 <- sin(theta1) * sin(theta2) * cos(theta3)
  m3 <- sin(theta1) * sin(theta2) * sin(theta3)

  cos_term <- x0 * m0 + x1 * m1 + x2 * m2 + x3 * m3
  log_raw <- cdm_log_density_core(
    log_fpt = log_fpt,
    tt = tt,
    a = a,
    v = v,
    sigma = sigma,
    sv = sv,
    cos_term = cos_term,
    dim = 4,
    log_surface = log(2 * pi^2)
  )

  dens <- exp(log_raw)
  dens[!is.finite(dens)] <- 0
  dens
}

# Joint density of (RT, R) under the Projected Spherical Diffusion Model
dPSDM <- function(rt, R, pars) {
  n <- length(rt)
  rt <- as.numeric(rt)
  R <- as.numeric(R)

  t0 <- pars[, "t0"]
  a <- pars[, "a"]
  v <- pars[, "v"]
  theta1 <- pars[, "theta1"] * pi
  sigma <- pars[, "s"]
  sv <- if ("sv" %in% colnames(pars)) as.numeric(pars[, "sv"]) else rep(0, n)

  tt <- pmax(rt - t0, 0)
  sig2 <- sigma * sigma
  v2 <- v * v
  v_scale <- sig2 / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02

  fpt_series <- series_sdm_fpt(tt, a = a, sigma = sigma)
  fpt_small <- v_scale * small_t_fpt_sdm(s, 1e-8 / (a * a))
  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  fpt <- (1 - w) * fpt_small + w * fpt_series
  fpt[!is.finite(fpt) | (fpt < 0)] <- 0
  log_fpt <- log(pmax(fpt, .Machine$double.xmin))

  c1 <- cos(theta1)
  s1a <- sin(theta1)
  cR <- cos(R)
  sR <- sin(R)
  A <- c1 * cR
  B <- s1a * sR

  log_raw <- rep(-Inf, n)
  idx0 <- (sv <= 0)
  if (any(idx0)) {
    k0 <- (a[idx0] * v[idx0]) / sig2[idx0]
    log_base0 <- log_fpt[idx0] - 0.5 * (v2[idx0] * tt[idx0]) / sig2[idx0]
    log_raw[idx0] <- log_base0 + k0 * A[idx0] + log_i0_stable(k0 * B[idx0]) - log(2)
  }

  idxv <- (sv > 0)
  if (any(idxv)) {
    ttv <- pmax(tt[idxv], .Machine$double.eps)
    sv2 <- sv[idxv]^2
    D <- sig2[idxv] + sv2 * ttv
    kv <- (a[idxv] * v[idxv]) / D
    log_basev <- log_fpt[idxv] +
      0.5 * 3 * (log(sig2[idxv]) - log(D)) +
      0.5 * (a[idxv]^2 * sv2) / (sig2[idxv] * D) -
      0.5 * (v2[idxv] * ttv) / D
    log_raw[idxv] <- log_basev + kv * A[idxv] + log_i0_stable(kv * B[idxv]) - log(2)
  }

  dens <- exp(log_raw)
  dens[!is.finite(dens)] <- 0
  dens
}

# Joint density of (RT, R, R2) under the Projected Hyper-Spherical Diffusion Model
dPHSDM <- function(rt, R, R2, pars) {
  n <- length(rt)
  rt <- as.numeric(rt)
  R <- as.numeric(R)
  R2 <- as.numeric(R2)

  t0 <- pars[, "t0"]
  a <- pars[, "a"]
  v <- pars[, "v"]
  theta1 <- pars[, "theta1"] * pi
  theta2 <- pars[, "theta2"] * pi
  sigma <- pars[, "s"]
  sv <- if ("sv" %in% colnames(pars)) as.numeric(pars[, "sv"]) else rep(0, n)

  tt <- pmax(rt - t0, 0)
  sig2 <- sigma * sigma
  v2 <- v * v
  v_scale <- sig2 / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02

  bessel <- load_bessel_values_hsdm()
  fpt_series <- series_hsdm_fpt(tt, a = a, sigma = sigma, bessel = bessel)
  fpt_small <- v_scale * small_t_fpt_hsdm(s, 1e-8 / (a * a), bessel = bessel)
  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  fpt <- (1 - w) * fpt_small + w * fpt_series
  fpt[!is.finite(fpt) | (fpt < 0)] <- 0
  log_fpt <- log(pmax(fpt, .Machine$double.xmin))

  st1 <- sin(theta1)
  ct1 <- cos(theta1)
  st2 <- sin(theta2)
  ct2 <- cos(theta2)
  sR <- sin(R)
  cR <- cos(R)
  sR2 <- sin(R2)
  cR2 <- cos(R2)

  A <- ct1 * cR + st1 * ct2 * sR * cR2
  B <- st1 * st2 * sR * sR2

  log_raw <- rep(-Inf, n)
  idx0 <- (sv <= 0)
  if (any(idx0)) {
    k0 <- (a[idx0] * v[idx0]) / sig2[idx0]
    log_base0 <- log_fpt[idx0] - 0.5 * (v2[idx0] * tt[idx0]) / sig2[idx0]
    log_raw[idx0] <- log_base0 + k0 * A[idx0] + log_i0_stable(k0 * B[idx0]) - log(pi)
  }

  idxv <- (sv > 0)
  if (any(idxv)) {
    ttv <- pmax(tt[idxv], .Machine$double.eps)
    sv2 <- sv[idxv]^2
    D <- sig2[idxv] + sv2 * ttv
    kv <- (a[idxv] * v[idxv]) / D
    log_basev <- log_fpt[idxv] +
      0.5 * 4 * (log(sig2[idxv]) - log(D)) +
      0.5 * (a[idxv]^2 * sv2) / (sig2[idxv] * D) -
      0.5 * (v2[idxv] * ttv) / D
    log_raw[idxv] <- log_basev + kv * A[idxv] + log_i0_stable(kv * B[idxv]) - log(pi)
  }

  dens <- exp(log_raw)
  dens[!is.finite(dens)] <- 0
  dens
}

log_likelihood_cdm_multi <- function(pars, dadm, model, min_ll = log(1e-10)) {
  like <- numeric(nrow(dadm))
  ok <- attr(pars, "ok")
  if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
  ok <- ok & (dadm$rt > pars[, "t0"])

  if (any(ok)) {
    if ("R3" %in% names(dadm)) {
      like[ok] <- dHSDM(dadm$rt[ok], dadm$R[ok], dadm$R2[ok], dadm$R3[ok], pars[ok, , drop = FALSE])
    } else if ("R2" %in% names(dadm)) {
      like[ok] <- dSDM(dadm$rt[ok], dadm$R[ok], dadm$R2[ok], pars[ok, , drop = FALSE])
    } else {
      like[ok] <- dCDM(dadm$rt[ok], dadm$R[ok], pars[ok, , drop = FALSE])
    }
  }
  like[is.na(like)] <- 0
  expand <- attr(dadm, "expand")
  if (is.null(expand)) expand <- seq_along(like)
  sum(pmax(min_ll, log(like[expand])))
}

log_likelihood_psdm <- function(pars, dadm, model, min_ll = log(1e-10)) {
  like <- numeric(nrow(dadm))
  ok <- attr(pars, "ok")
  if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
  ok <- ok & (dadm$rt > pars[, "t0"])
  if (any(ok)) {
    like[ok] <- dPSDM(dadm$rt[ok], dadm$R[ok], pars[ok, , drop = FALSE])
  }
  like[is.na(like)] <- 0
  expand <- attr(dadm, "expand")
  if (is.null(expand)) expand <- seq_along(like)
  sum(pmax(min_ll, log(like[expand])))
}

log_likelihood_phsdm <- function(pars, dadm, model, min_ll = log(1e-10)) {
  like <- numeric(nrow(dadm))
  ok <- attr(pars, "ok")
  if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
  ok <- ok & (dadm$rt > pars[, "t0"])
  if (any(ok)) {
    like[ok] <- dPHSDM(dadm$rt[ok], dadm$R[ok], dadm$R2[ok], pars[ok, , drop = FALSE])
  }
  like[is.na(like)] <- 0
  expand <- attr(dadm, "expand")
  if (is.null(expand)) expand <- seq_along(like)
  sum(pmax(min_ll, log(like[expand])))
}

#' The Spherical Diffusion Model (SDM)
#'
#' Model file to estimate the Spherical Diffusion Model (SDM) in EMC2.
#'
#' @details
#' Parameters and defaults are on the transformed scale.
#' All parameters not specified in the `formula` of `design()` are
#' assumed constant with the defaults in `p_types`.
#'
#' | Parameter | Transform | Natural scale | Default   | Interpretation |
#' |-----------|-----------|---------------|-----------|----------------|
#' | v         | log       | [0, Inf)      | log(1)    | Drift magnitude |
#' | theta1    | pnorm     | [0, 1)        | qnorm(.5) | Drift polar angle, mapped to 0..pi |
#' | theta2    | pnorm     | [0, 1)        | qnorm(.5) | Drift azimuth angle, mapped to -pi..pi |
#' | a         | log       | [0, Inf)      | log(1)    | Boundary radius |
#' | t0        | log       | [0, Inf)      | log(0)    | Non-decision time |
#' | s         | log       | [0, Inf)      | log(1)    | Diffusion scale |
#' | sv        | log       | [0, Inf)      | log(0)    | Drift SD across trials |
#'
#' The response variables are `rt`, `R`, and `R2`.
#' This implementation assumes fixed boundaries (no decay/collapsing boundary parameter).
#'
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @export
SDM <- function() {
  list(
    type = "CDM",
    c_name = "CDM",
    p_types = c(
      "v" = log(1),
      "theta1" = qnorm(0.5),
      "theta2" = qnorm(0.5),
      "a" = log(1),
      "t0" = log(0),
      "s" = log(1),
      "sv" = log(0)
    ),
    transform = list(func = c(
      v = "exp",
      theta1 = "pnorm",
      theta2 = "pnorm",
      a = "exp",
      t0 = "exp",
      s = "exp",
      sv = "exp"
    )),
    bound = list(
      minmax = cbind(
        v = c(0, 10),
        theta1 = c(1e-4, 1 - 1e-4),
        theta2 = c(1e-4, 1 - 1e-4),
        a = c(1e-4, 7),
        t0 = c(0.05, Inf),
        s = c(1e-4, Inf),
        sv = c(0, 2)
      ),
      exception = c(sv = 0, t0 = 0)
    ),
    Ttransform = function(pars, dadm) { pars },
    rfun = function(data = NULL, pars) rSDM(pars, attr(pars, "ok")),
    dfun = function(rt, R, R2, pars) dSDM(rt, R, R2, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_cdm_multi(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

#' The Hyper-Spherical Diffusion Model (HSDM)
#'
#' Model file to estimate the Hyper-Spherical Diffusion Model (HSDM) in EMC2.
#'
#' @details
#' Parameters and defaults are on the transformed scale.
#' All parameters not specified in the `formula` of `design()` are
#' assumed constant with the defaults in `p_types`.
#'
#' | Parameter | Transform | Natural scale | Default   | Interpretation |
#' |-----------|-----------|---------------|-----------|----------------|
#' | v         | log       | [0, Inf)      | log(1)    | Drift magnitude |
#' | theta1    | pnorm     | [0, 1)        | qnorm(.5) | Drift angle, mapped to 0..pi |
#' | theta2    | pnorm     | [0, 1)        | qnorm(.5) | Drift angle, mapped to 0..pi |
#' | theta3    | pnorm     | [0, 1)        | qnorm(.5) | Drift angle, mapped to -pi..pi |
#' | a         | log       | [0, Inf)      | log(1)    | Boundary radius |
#' | t0        | log       | [0, Inf)      | log(0)    | Non-decision time |
#' | s         | log       | [0, Inf)      | log(1)    | Diffusion scale |
#' | sv        | log       | [0, Inf)      | log(0)    | Drift SD across trials |
#'
#' The response variables are `rt`, `R`, `R2`, and `R3`.
#' This implementation assumes fixed boundaries (no decay/collapsing boundary parameter).
#'
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @export
HSDM <- function() {
  list(
    type = "CDM",
    c_name = "CDM",
    p_types = c(
      "v" = log(1),
      "theta1" = qnorm(0.5),
      "theta2" = qnorm(0.5),
      "theta3" = qnorm(0.5),
      "a" = log(1),
      "t0" = log(0),
      "s" = log(1),
      "sv" = log(0)
    ),
    transform = list(func = c(
      v = "exp",
      theta1 = "pnorm",
      theta2 = "pnorm",
      theta3 = "pnorm",
      a = "exp",
      t0 = "exp",
      s = "exp",
      sv = "exp"
    )),
    bound = list(
      minmax = cbind(
        v = c(0, 10),
        theta1 = c(1e-4, 1 - 1e-4),
        theta2 = c(1e-4, 1 - 1e-4),
        theta3 = c(1e-4, 1 - 1e-4),
        a = c(1e-4, 7),
        t0 = c(0.05, Inf),
        s = c(1e-4, Inf),
        sv = c(0, 2)
      ),
      exception = c(sv = 0, t0 = 0)
    ),
    Ttransform = function(pars, dadm) { pars },
    rfun = function(data = NULL, pars) rHSDM(pars, attr(pars, "ok")),
    dfun = function(rt, R, R2, R3, pars) dHSDM(rt, R, R2, R3, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_cdm_multi(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

#' The Projected Spherical Diffusion Model (PSDM)
#'
#' Model file to estimate the Projected Spherical Diffusion Model (PSDM) in EMC2.
#'
#' @details
#' Parameters and defaults are on the transformed scale.
#' All parameters not specified in the `formula` of `design()` are
#' assumed constant with the defaults in `p_types`.
#'
#' | Parameter | Transform | Natural scale | Default   | Interpretation |
#' |-----------|-----------|---------------|-----------|----------------|
#' | v         | log       | [0, Inf)      | log(1)    | Drift magnitude |
#' | theta1    | pnorm     | [0, 1)        | qnorm(.5) | Drift polar angle, mapped to 0..pi |
#' | a         | log       | [0, Inf)      | log(1)    | Boundary radius |
#' | t0        | log       | [0, Inf)      | log(0)    | Non-decision time |
#' | s         | log       | [0, Inf)      | log(1)    | Diffusion scale |
#' | sv        | log       | [0, Inf)      | log(0)    | Drift SD across trials |
#'
#' The response variables are `rt` and `R`.
#' This implementation assumes fixed boundaries (no decay/collapsing boundary parameter).
#'
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @export
PSDM <- function() {
  list(
    type = "CDM",
    c_name = "PSDM",
    p_types = c(
      "v" = log(1),
      "theta1" = qnorm(0.5),
      "a" = log(1),
      "t0" = log(0),
      "s" = log(1),
      "sv" = log(0)
    ),
    transform = list(func = c(
      v = "exp",
      theta1 = "pnorm",
      a = "exp",
      t0 = "exp",
      s = "exp",
      sv = "exp"
    )),
    bound = list(
      minmax = cbind(
        v = c(0, 10),
        theta1 = c(1e-4, 1 - 1e-4),
        a = c(1e-4, 7),
        t0 = c(0.05, Inf),
        s = c(1e-4, Inf),
        sv = c(0, 2)
      ),
      exception = c(sv = 0, t0 = 0)
    ),
    Ttransform = function(pars, dadm) { pars },
    rfun = function(data = NULL, pars) rPSDM(pars, attr(pars, "ok")),
    dfun = function(rt, R, pars) dPSDM(rt, R, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_psdm(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

#' The Projected Hyper-Spherical Diffusion Model (PHSDM)
#'
#' Model file to estimate the Projected Hyper-Spherical Diffusion Model (PHSDM) in EMC2.
#'
#' @details
#' Parameters and defaults are on the transformed scale.
#' All parameters not specified in the `formula` of `design()` are
#' assumed constant with the defaults in `p_types`.
#'
#' | Parameter | Transform | Natural scale | Default   | Interpretation |
#' |-----------|-----------|---------------|-----------|----------------|
#' | v         | log       | [0, Inf)      | log(1)    | Drift magnitude |
#' | theta1    | pnorm     | [0, 1)        | qnorm(.5) | Drift angle, mapped to 0..pi |
#' | theta2    | pnorm     | [0, 1)        | qnorm(.5) | Drift angle, mapped to 0..pi |
#' | a         | log       | [0, Inf)      | log(1)    | Boundary radius |
#' | t0        | log       | [0, Inf)      | log(0)    | Non-decision time |
#' | s         | log       | [0, Inf)      | log(1)    | Diffusion scale |
#' | sv        | log       | [0, Inf)      | log(0)    | Drift SD across trials |
#'
#' The response variables are `rt`, `R`, and `R2`.
#' This implementation assumes fixed boundaries (no decay/collapsing boundary parameter).
#'
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @export
PHSDM <- function() {
  list(
    type = "CDM",
    c_name = "PHSDM",
    p_types = c(
      "v" = log(1),
      "theta1" = qnorm(0.5),
      "theta2" = qnorm(0.5),
      "a" = log(1),
      "t0" = log(0),
      "s" = log(1),
      "sv" = log(0)
    ),
    transform = list(func = c(
      v = "exp",
      theta1 = "pnorm",
      theta2 = "pnorm",
      a = "exp",
      t0 = "exp",
      s = "exp",
      sv = "exp"
    )),
    bound = list(
      minmax = cbind(
        v = c(0, 10),
        theta1 = c(1e-4, 1 - 1e-4),
        theta2 = c(1e-4, 1 - 1e-4),
        a = c(1e-4, 7),
        t0 = c(0.05, Inf),
        s = c(1e-4, Inf),
        sv = c(0, 2)
      ),
      exception = c(sv = 0, t0 = 0)
    ),
    Ttransform = function(pars, dadm) { pars },
    rfun = function(data = NULL, pars) rPHSDM(pars, attr(pars, "ok")),
    dfun = function(rt, R, R2, pars) dPHSDM(rt, R, R2, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_phsdm(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
