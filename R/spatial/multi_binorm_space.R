# Script fits range of spatial models using the bi-normal curve as the basis for
# multifilament catch data. Candidate models include spatially varying
# deviations for all parameters, only means, only variance, and similar
# combinations where means/variance parameters spatially vary, but univariate
# random effects assumed for other parameters.

library(RTMB)

# Read in multifilament catch data
data <- read.csv("data/multi_catch_net_utm.csv")

# clean up data a bit
names(data)[17] <- "easting_km"
names(data)[18] <- "northing_km"

n_sites <- length(unique(data$net_id))
catches <- as.matrix(data[, which(grepl("mesh", colnames(data)))])
colnames(catches) <- NULL

gn_meshes <- seq(2.0, 5.0, by = 0.25) # range of stretch measure meshes in inches
gn_meshes <- gn_meshes * 25.4 # convert to mm

locs <- unique(data[, c("easting_km", "northing_km")])
plot(locs)

# create mesh to model spatial correlation among sites
mesh <- INLA::inla.mesh.2d(loc = locs, max.edge = c(6, 50))

plot(mesh)
points(locs, col = "#075057", pch = 16)
mesh$n
n_sites

data <- within(data, net <- as.numeric(interaction(data$net_id, drop = TRUE, lex.order = F)))
data <- data[order(data$net), ]

# Create list for RTMB
data <- list(
  meshes = gn_meshes,
  n_site = length(unique(data$net)),
  n_obs = nrow(data),
  loc_idx = data$net,
  lens = data[, 2],
  catches = catches
)
data$rel_size <- data$meshes / data$meshes[1]

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)

data$spde <- spde$param.inla[c("M0", "M1", "M2")]

#### all parameters spatially varying ####

# Fit once and re-fit using estimates below

pars <- list(
  ln_k1 = 5.5,
  ln_k2 = 3.1,
  ln_k3 = 5.8,
  ln_k4 = 4.1,
  ln_c = -0.1,
  eps_k1 = rep(0, mesh$n),
  eps_k2 = rep(0, mesh$n),
  eps_k3 = rep(0, mesh$n),
  eps_k4 = rep(0, mesh$n),
  eps_c = rep(0, mesh$n),
  log_tau1 = -2,
  log_tau2 = 1,
  log_tau3 = -2,
  log_tau4 = -1,
  log_tau5 = 1,
  log_kappa1 = -1,
  log_kappa2 = -1,
  log_kappa3 = -1,
  log_kappa4 = -1,
  log_kappa5 = -1
)

# Update starting values
pars <- list(
  ln_k1 = rep$par.fixed[1],
  ln_k2 = rep$par.fixed[2],
  ln_k3 = rep$par.fixed[3],
  ln_k4 = rep$par.fixed[4],
  ln_c = rep$par.fixed[5],
  eps_k1 = rep$par.random[1:mesh$n],
  eps_k2 = rep$par.random[(mesh$n + 1):(2 * mesh$n)],
  eps_k3 = rep(0, mesh$n), # rep$par.random[(2*mesh$n + 1):(3*mesh$n)],
  eps_k4 = rep$par.random[(3 * mesh$n + 1):(4 * mesh$n)],
  eps_c = rep$par.random[(4 * mesh$n + 1):(5 * mesh$n)],
  log_tau1 = rep$par.fixed[6],
  log_tau2 = rep$par.fixed[7],
  log_tau3 = 1.5, # rep$par.fixed[8],
  log_tau4 = rep$par.fixed[9],
  log_tau5 = rep$par.fixed[10],
  log_kappa1 = rep$par.fixed[11],
  log_kappa2 = rep$par.fixed[12],
  log_kappa3 = -.1, # rep$par.fixed[13],
  log_kappa4 = rep$par.fixed[14],
  log_kappa5 = rep$par.fixed[15]
)

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

f <- function(parameters) {
  getAll(data, parameters)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  tau3 <- exp(log_tau3)
  tau4 <- exp(log_tau4)
  tau5 <- exp(log_tau5)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)
  kappa3 <- exp(log_kappa3)
  kappa4 <- exp(log_kappa4)
  kappa5 <- exp(log_kappa5)

  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)
  Q3 <- Q_spde(spde, kappa3)
  Q4 <- Q_spde(spde, kappa4)
  Q5 <- Q_spde(spde, kappa5)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k1, 0, Q1, log = TRUE, scale = 1 / tau1)
  jnll <- jnll - dgmrf(eps_k2, 0, Q2, log = TRUE, scale = 1 / tau2)
  jnll <- jnll - dgmrf(eps_k3, 0, Q3, log = TRUE, scale = 1 / tau3)
  jnll <- jnll - dgmrf(eps_k4, 0, Q4, log = TRUE, scale = 1 / tau4)
  jnll <- jnll - dgmrf(eps_c, 0, Q5, log = TRUE, scale = 1 / tau5)

  k1 <- exp(ln_k1 + eps_k1)
  k2 <- exp(ln_k2 + eps_k2)
  k3 <- exp(ln_k3 + eps_k3)
  k4 <- exp(ln_k4 + eps_k4)
  c <- exp(ln_c + eps_c)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[loc_idx[i]] * rel_size[j])^2 /
        (2 * k2[loc_idx[i]]^2 * rel_size[j]^2)) +
        (c[loc_idx[i]] * exp(-(lens[i] - k3[loc_idx[i]] * rel_size[j])^2 /
          (2 * k4[loc_idx[i]]^2 * rel_size[j]^2)))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))

  # Derived variables

  range1 <- sqrt(8) / exp(log_kappa1)
  range2 <- sqrt(8) / exp(log_kappa2)
  range3 <- sqrt(8) / exp(log_kappa3)
  range4 <- sqrt(8) / exp(log_kappa4)
  range5 <- sqrt(8) / exp(log_kappa5)

  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k2 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))
  sig_o_k3 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa3))
  sig_o_k4 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa4))
  sig_o_k5 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa5))

  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(range3)
  ADREPORT(range4)
  ADREPORT(range5)
  ADREPORT(sig_o_k1)
  ADREPORT(sig_o_k2)
  ADREPORT(sig_o_k3)
  ADREPORT(sig_o_k4)
  ADREPORT(sig_o_k5)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_k1", "eps_k2", "eps_k3", "eps_k4", "eps_c"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

n_newton <- 3 # number of Newton steps
tryCatch(
  for (n in 1:n_newton) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite optults
    opt <- nlminb(new_par, obj$fn, obj$gr)
  },
  error = function(e) {
    err <<- conditionMessage(e)
  }
)
rep <- sdreport(obj)
rep

# Large SEs, but meets convergence criteria
# Fit with TMBhelper

TMBhelper::fit_tmb(obj = obj, newtonsteps = 3, getReportCovariance = TRUE)

#### Means spatially varying, others fixed ####
pars <- list(
  ln_k1 = 5.5,
  ln_k2 = 3.1,
  ln_k3 = 5.8,
  ln_k4 = 4.1,
  ln_c = -0.1,
  eps_k1 = rep(0, mesh$n),
  eps_k3 = rep(0, mesh$n),
  log_tau1 = -2, # for k1
  log_tau2 = 1, # for k3
  log_kappa1 = -1, # for k1
  log_kappa2 = -1 # for k3
)

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

f <- function(parameters) {
  getAll(data, parameters)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)

  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k1, 0, Q1, log = TRUE, scale = 1 / tau1)
  jnll <- jnll - dgmrf(eps_k3, 0, Q2, log = TRUE, scale = 1 / tau2)


  k1 <- exp(ln_k1 + eps_k1)
  k2 <- exp(ln_k2)
  k3 <- exp(ln_k3 + eps_k3)
  k4 <- exp(ln_k4)
  c <- exp(ln_c)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[loc_idx[i]] * rel_size[j])^2 /
        (2 * k2^2 * rel_size[j]^2)) +
        (c * exp(-(lens[i] - k3[loc_idx[i]] * rel_size[j])^2 /
          (2 * k4^2 * rel_size[j]^2)))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))

  # Derived variables
  range1 <- sqrt(8) / exp(log_kappa1)
  range2 <- sqrt(8) / exp(log_kappa2)

  #
  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k3 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  #
  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_k1)
  ADREPORT(sig_o_k3)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_k1", "eps_k3"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

n_newton <- 3 # number of Newton steps
tryCatch(
  for (n in 1:n_newton) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite optults
    opt <- nlminb(new_par, obj$fn, obj$gr)
  },
  error = function(e) {
    err <<- conditionMessage(e)
  }
)
rep <- sdreport(obj)
rep

# No issues with large SEs
fit_means <- TMBhelper::fit_tmb(obj = obj, newtonsteps = 3, getReportCovariance = TRUE)

#### Variances spatial, everything else fixed ####

pars <- list(
  ln_k1 = 5.5,
  ln_k2 = 3.1,
  ln_k3 = 5.8,
  ln_k4 = 4.1,
  ln_c = -0.1,
  eps_k2 = rep(0, mesh$n),
  eps_k4 = rep(0, mesh$n),
  log_tau1 = -2, # for k2
  log_tau2 = 1, # for k4
  log_kappa1 = -1, # for k2
  log_kappa2 = -1 # for k4
)

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

f <- function(parameters) {
  getAll(data, parameters)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)

  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k2, 0, Q1, log = TRUE, scale = 1 / tau1)
  jnll <- jnll - dgmrf(eps_k4, 0, Q2, log = TRUE, scale = 1 / tau2)

  k1 <- exp(ln_k1)
  k2 <- exp(ln_k2 + eps_k2)
  k3 <- exp(ln_k3)
  k4 <- exp(ln_k4 + eps_k4)
  c <- exp(ln_c)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1 * rel_size[j])^2 /
        (2 * k2[loc_idx[i]]^2 * rel_size[j]^2)) +
        (c * exp(-(lens[i] - k3 * rel_size[j])^2 /
          (2 * k4[loc_idx[i]]^2 * rel_size[j]^2)))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))

  # Derived variables
  range1 <- sqrt(8) / exp(log_kappa1) # k2 = variance 1st curve
  range2 <- sqrt(8) / exp(log_kappa2) # k4 = variance 2nd curve

  #
  sig_o_k2 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k4 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  #
  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_k2)
  ADREPORT(sig_o_k4)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_k2", "eps_k4"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

n_newton <- 3 # number of Newton steps
tryCatch(
  for (n in 1:n_newton) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite optults
    opt <- nlminb(new_par, obj$fn, obj$gr)
  },
  error = function(e) {
    err <<- conditionMessage(e)
  }
)
rep <- sdreport(obj)
rep

# NA/NaN function evaluation
TMBhelper::fit_tmb(obj = obj, newtonsteps = 3, getReportCovariance = TRUE)

#### Variances spatially varying, IID everything else ####
pars <- list(
  ln_k1 = 5.5,
  ln_k2 = 3.1,
  ln_k3 = 5.8,
  ln_k4 = 4.1,
  ln_c = -0.1,
  eps_k1 = rep(0, n_sites),
  eps_k2 = rep(0, mesh$n),
  eps_k3 = rep(0, n_sites),
  eps_k4 = rep(0, mesh$n),
  eps_c = rep(0, n_sites),
  log_tau1 = -2, # for k2
  log_tau2 = 1, # for k4
  log_kappa1 = -1, # for k2
  log_kappa2 = -1, # for k4
  ln_sd_k1 = 0.1,
  ln_sd_k3 = 0.1,
  ln_sd_c = 0.1
)

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

f <- function(parameters) {
  getAll(data, parameters)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)

  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k2, 0, Q1, log = TRUE, scale = 1 / tau1) # mean0 k1 deviates
  jnll <- jnll - sum(dnorm(eps_k1, 0, exp(ln_sd_k1), TRUE)) # k2
  jnll <- jnll - dgmrf(eps_k4, 0, Q2, log = TRUE, scale = 1 / tau2) # mean0 k3 deviates
  jnll <- jnll - sum(dnorm(eps_k3, 0, exp(ln_sd_k3), TRUE)) # k4
  jnll <- jnll - sum(dnorm(eps_c, 0, exp(ln_sd_c), TRUE)) # c

  k1 <- exp(ln_k1 + eps_k1) # spatially correlated
  k2 <- exp(ln_k2 + eps_k2) # IID
  k3 <- exp(ln_k3 + eps_k3) # spatially correlated
  k4 <- exp(ln_k4 + eps_k4) # IID
  c <- exp(ln_c + eps_c) # IID

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[loc_idx[i]] * rel_size[j])^2 /
        (2 * k2[loc_idx[i]]^2 * rel_size[j]^2)) +
        (c[loc_idx[i]] * exp(-(lens[i] - k3[loc_idx[i]] * rel_size[j])^2 /
          (2 * k4[loc_idx[i]]^2 * rel_size[j]^2)))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))

  # Derived variables
  range1 <- sqrt(8) / exp(log_kappa1)
  range2 <- sqrt(8) / exp(log_kappa2)

  #
  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k3 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  #
  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_k1)
  ADREPORT(sig_o_k3)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_k1", "eps_k2", "eps_k3", "eps_k4", "eps_c"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

n_newton <- 3 # number of Newton steps
tryCatch(
  for (n in 1:n_newton) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite optults
    opt <- nlminb(new_par, obj$fn, obj$gr)
  },
  error = function(e) {
    err <<- conditionMessage(e)
  }
)
rep <- sdreport(obj)
rep

# Poorly estimates spatial field
fit_vars <- TMBhelper::fit_tmb(obj = obj, newtonsteps = 3, getReportCovariance = TRUE)

#### Means spatially varying, IID variances and c ####

pars <- list(
  ln_k1 = 5.5,
  ln_k2 = 3.1,
  ln_k3 = 5.8,
  ln_k4 = 4.1,
  ln_c = -0.1,
  eps_k1 = rep(0, mesh$n),
  eps_k2 = rep(0, n_sites),
  eps_k3 = rep(0, mesh$n),
  eps_k4 = rep(0, n_sites),
  eps_c = rep(0, n_sites),
  log_tau1 = -2, # for k1
  log_tau2 = 1, # for k3
  log_kappa1 = -1, # for k1
  log_kappa2 = -1, # for k3
  ln_sd_k2 = 0.1,
  ln_sd_k4 = 0.1,
  ln_sd_c = 0.1
)

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

f <- function(parameters) {
  getAll(data, parameters)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)

  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k1, 0, Q1, log = TRUE, scale = 1 / tau1) # mean0 k1 deviates
  jnll <- jnll - sum(dnorm(eps_k2, 0, exp(ln_sd_k2), TRUE)) # k2
  jnll <- jnll - dgmrf(eps_k3, 0, Q2, log = TRUE, scale = 1 / tau2) # mean0 k3 deviates
  jnll <- jnll - sum(dnorm(eps_k4, 0, exp(ln_sd_k4), TRUE)) # k4
  jnll <- jnll - sum(dnorm(eps_c, 0, exp(ln_sd_c), TRUE)) # c

  k1 <- exp(ln_k1 + eps_k1) # spatially correlated
  k2 <- exp(ln_k2 + eps_k2) # IID
  k3 <- exp(ln_k3 + eps_k3) # spatially correlated
  k4 <- exp(ln_k4 + eps_k4) # IID
  c <- exp(ln_c + eps_c) # IID

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[loc_idx[i]] * rel_size[j])^2 /
        (2 * k2[loc_idx[i]]^2 * rel_size[j]^2)) +
        (c[loc_idx[i]] * exp(-(lens[i] - k3[loc_idx[i]] * rel_size[j])^2 /
          (2 * k4[loc_idx[i]]^2 * rel_size[j]^2)))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))

  # Derived variables
  range1 <- sqrt(8) / exp(log_kappa1)
  range2 <- sqrt(8) / exp(log_kappa2)

  #
  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k3 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  #
  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_k1)
  ADREPORT(sig_o_k3)

  jnll
}


obj <- MakeADFun(f, pars, random = c("eps_k1", "eps_k2", "eps_k3", "eps_k4", "eps_c"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

n_newton <- 3 # number of Newton steps
tryCatch(
  for (n in 1:n_newton) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite optults
    opt <- nlminb(new_par, obj$fn, obj$gr)
  },
  error = function(e) {
    err <<- conditionMessage(e)
  }
)
rep <- sdreport(obj)
rep

fit_means <- TMBhelper::fit_tmb(obj = obj, newtonsteps = 3, getReportCovariance = TRUE)
fit_means$objective # 17407.66
fit_means$AIC # 34839.33

#### End of script ####
