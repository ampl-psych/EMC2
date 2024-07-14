rm(list = ls())
devtools::load_all()

N <- 10

# loading matrix
Lambda <- matrix(c(1, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1,
                   .7, 0, 0, 0,
                   0, .5, 0, 0,
                   0, 0, .6, 0,
                   0, 0, 0, .5,
                   .9, 0, 0, 0,
                   0, .8, 0, 0,
                   0, 0, .4, 0,
                   0, 0, 0, .3
), nrow = 12, ncol =4, byrow = T)

K <- matrix(c(0, 0,
              0, 0,
              0, 0,
              0, 0,
              0, 0,
              0, 0,
              0, 0,
              0, 0,
              0.5, 0,
              0.3, 0,
              0.6, 0,
              0, 0), ncol = 2, byrow = T)

# K <- matrix(0, nrow = 12, ncol = 2)


G <- matrix(c(0, 0.3,
              0, 0,
              0, 0,
              0, 0), ncol = 2, byrow = T)
# G <- matrix(0, nrow = 3, ncol = 2)

B <- matrix(c(0, -.5, .5, 0,
              0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0), nrow = 4, ncol = 4, byrow = T)

# Some fixed covariate
x <- mvtnorm::rmvnorm(N, mean = c(0, 0), sigma = diag(2)/10)
colnames(x) <- paste0("v", 1:(ncol(x)))

lat_cov <- matrix(c(.3, 0, 0, 0,
                    0, .2, 0, 0.1,
                    0, 0, .1, 0.05,
                    0, 0.1, 0.05, .2), nrow = 4)
# loadings with no latent predictors on them
xi <- x %*% t(G[-1,]) + mvtnorm::rmvnorm(N, mean = c(0, 0, 0), sigma = lat_cov[2:4, 2:4])

zeta <- rnorm(N, sd = sqrt(lat_cov[1,1]))
eta1 <- x %*% t(G[1,, drop = F]) + tcrossprod(xi, t(B[1,2:4])) + zeta
eta <- cbind(eta1, xi)

# error term
#
# debug(sampled_p_vector)
design <- design(
  factors =list(subjects=1:N,S=c("right", "left")),
  Rlevels=c("right", "left"),matchfun=function(d)d$S==d$lR,
  formula=list(m~lM,s~1,t0~1),
  model=LNR)


Theta <- diag(0.1, nrow = 12)
epsilon <- mvtnorm::rmvnorm(N, mean = rep(0, ncol(Theta)), sigma = Theta)
mu <- rep(c(.5, -.3, log(.4), log(.2)), 3)
dat <- matrix(mu, nrow = N, ncol = 12, byrow = T) + x %*% t(K) + tcrossprod(eta, Lambda) + epsilon

save(dat, file = "random_effects.RData")

colnames(dat) <- rep(names(sampled_p_vector(design)), times = 3)
rownames(dat) <- as.character(1:N)

all_data <- list()
all_data[[1]] <- make_data(parameters =  as.matrix(dat[,1:4]),
                           design = design, n_trials = 100)
all_data[[2]] <- make_data(parameters =  as.matrix(dat[,5:8]),
                           design = design, n_trials = 100)
all_data[[3]] <- make_data(parameters =  as.matrix(dat[,9:12]),
                           design = design, n_trials = 100)

Lambda_mat <- as.matrix(Lambda)
Lambda_mat[Lambda_mat != 1 & Lambda_mat != 0] <- Inf #Inf means that the parameter is freely estimated
B_mat <- as.matrix(B)
B_mat[B_mat != 0] <- Inf

K_mat <- as.matrix(K)
K_mat[K_mat != 0] <- Inf

G_mat <- as.matrix(G)
G_mat[G_mat != 0] <- Inf

prior_S <- prior(design = list(design, design, design), type = "SEM", Lambda_mat = Lambda_mat, B_mat = B_mat,
                 K_mat = K_mat, G_mat = G_mat, x = x)
plot_prior(prior_S, list(design, design, design), Lambda_mat = Lambda_mat, B_mat = B_mat,
           K_mat = K_mat, G_mat = G_mat, covariates = x,
           selection = "correlation", N = 1e4)

# load("SEM_LNR_with_G_small.RData")

# x <- samplers[[1]]$xeta

devtools::load_all()
samplers <- make_emc(data = all_data, design = rep(list(design), 3), type = "SEM",
                     Lambda_mat = Lambda_mat, B_mat = B_mat, K_mat = K_mat, G_mat = G_mat,
                     covariates= x)

fit(samplers, cores_per_chain = 4, fileName = "testSem.RData")


