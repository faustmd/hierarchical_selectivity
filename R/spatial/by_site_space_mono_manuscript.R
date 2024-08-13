# Fit spatial version of by-site hierarchical models to
# monofilament catch data from gear comparison study 2007-2008 & 2011-2013.
# Currently set up to allow all parameters to vary by site for each model.

library(RTMB)

# Read in monofilament catch data
data <- readRDS("data/by_net_mono.rds")

# remove this site because coordinates are unreliable
data <- data[data$CRNII != 2012005, ]

# convert coordinates from decimal degrees to UTM (km)
library(magrittr)

data %>%
  sf::st_as_sf(coords = c("LONG", "LAT")) %>%
  sf::st_set_crs(4326) %>%
  sf::st_transform(32617) %>%
  sf::st_coordinates() -> utm_coords

utm_coords_km <- utm_coords / 1000 # convert to km for ease of fitting

data <- cbind(data, utm_coords_km)

names(data) <- c(
  "net_id", "length_class", "mesh_1.5", "mesh_2.0", "mesh_2.5",
  "mesh_3.0", "mesh_3.5", "mesh_4.0", "mesh_4.5", "mesh_5.0",
  "mesh_5.5", "mesh_6.0", "mesh_6.5", "mesh_7.0", "latitude",
  "longitude", "easting_km", "northing_km"
)

n_sites <- length(unique(data$net_id))
catches <- as.matrix(data[, which(grepl("mesh", colnames(data)))])
colnames(catches) <- NULL

locs <- unique(data[, c("easting_km", "northing_km")])
plot(locs)

# create mesh to model spatial correlation among sites
mesh <- INLA::inla.mesh.2d(loc = locs, max.edge = c(6, 50))

plot(mesh)
points(locs, col = "#075057", pch = 16)
mesh$n # this is number of nodes (= random effects) in the mesh
n_sites # number of sites in data

data <- within(data, net <- as.numeric(interaction(data$net_id, drop = TRUE, lex.order = F)))
data <- data[order(data$net), ]

gn_meshes <- seq(1.5, 7.0, by = 0.5) # range of stretch measure meshes in inches
gn_meshes <- gn_meshes * 25.4 # convert to mm

# Create a list for RTMB
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

#### Normal location ####

pars <- list(
  ln_k = log(225),
  ln_sigma = log(75),
  eps_k = rep(0, mesh$n),
  eps_sigma = rep(0, mesh$n),
  log_tau1 = 2.12,
  log_tau2 = -2.0,
  log_kappa1 = -1.43,
  log_kappa2 = 0
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
  catches <- OBS(catches)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)

  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k, 0, Q1, log = TRUE, scale = 1 / tau1)
  jnll <- jnll - dgmrf(eps_sigma, 0, Q2, log = TRUE, scale = 1 / tau2)

  k <- exp(ln_k + eps_k)
  sigma <- exp(ln_sigma + eps_sigma)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k[loc_idx[i]] * rel_size[j])^2 /
        (2 * sigma[loc_idx[i]]^2))
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

  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k2 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  ADREPORT(range1)
  ADREPORT(range2)
  # ADREPORT(sig_o_k1)
  # ADREPORT(sig_o_k2)

  jnll
}

# Fit the model

obj <- MakeADFun(f, pars, random = c("eps_k", "eps_sigma"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

n_newton <- 3 # number of Newton steps
tryCatch(
  for (n in 1:n_newton) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite opt
    opt <- nlminb(new_par, obj$fn, obj$gr)
  },
  error = function(e) {
    err <<- conditionMessage(e)
  }
)

rep <- sdreport(obj)
rep

#### Normal ####

# starting parameter values
pars <- list(
  ln_k1 = 5.2,
  ln_k2 = 3.8,
  eps_k1 = rep(0, mesh$n),
  eps_k2 = rep(0, mesh$n),
  log_tau1 = 2.12,
  log_tau2 = -2.0,
  log_kappa1 = -1.43,
  log_kappa2 = 0
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
  catches <- OBS(catches)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)

  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k1, 0, Q1, log = TRUE, scale = 1 / tau1)
  jnll <- jnll - dgmrf(eps_k2, 0, Q2, log = TRUE, scale = 1 / tau2)

  k1 <- exp(ln_k1 + eps_k1)
  k2 <- exp(ln_k2 + eps_k2)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[loc_idx[i]] * rel_size[j])^2 /
        (2 * k2[loc_idx[i]]^2 * rel_size[j]^2))
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

  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k2 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_k1)
  ADREPORT(sig_o_k2)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_k1", "eps_k2"))

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$objective
n_newton <- 3 # number of Newton steps
tryCatch(
  for (n in 1:n_newton) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite opt
    opt <- nlminb(new_par, obj$fn, obj$gr)
  },
  error = function(e) {
    err <<- conditionMessage(e)
  }
)
rep <- sdreport(obj)
rep
#### Gamma ####

pars <- list(
  ln_alpha = 3,
  ln_k = 1.5,
  eps_alpha = rep(0, mesh$n),
  eps_k = rep(0, mesh$n),
  log_tau1 = -1, # 2.2,
  log_tau2 = -1,
  log_kappa1 = 1,
  log_kappa2 = 1
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
  jnll <- jnll - dgmrf(eps_alpha, 0, Q1, log = TRUE, scale = 1 / tau1)
  jnll <- jnll - dgmrf(eps_k, 0, Q2, log = TRUE, scale = 1 / tau2)

  alpha <- exp(ln_alpha + eps_alpha)
  k <- exp(ln_k + eps_k)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- ((lens[i] / ((alpha[loc_idx[i]] - 1) * k[loc_idx[i]] * rel_size[j]))^(alpha[loc_idx[i]] - 1)) *
        (exp(alpha[loc_idx[i]] - 1 - (lens[i] / (k[loc_idx[i]] * rel_size[j]))))
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

  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k2 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_k1)
  ADREPORT(sig_o_k2)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_alpha", "eps_k"))
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

#### Lognormal ####

pars <- list(
  ln_mu = 5.2,
  ln_sigma = 0.1,
  eps_mu = rep(0, mesh$n),
  eps_sigma = rep(0, mesh$n),
  log_tau1 = 2.5,
  log_tau2 = -.1,
  log_kappa1 = -1.13,
  log_kappa2 = 0.2
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
  jnll <- jnll - dgmrf(eps_mu, 0, Q1, log = TRUE, scale = 1 / tau1)
  jnll <- jnll - dgmrf(eps_sigma, 0, Q2, log = TRUE, scale = 1 / tau2)

  mu <- (ln_mu + eps_mu)
  sigma <- (ln_sigma + eps_sigma)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- (rel_size[j] / lens[i]) *
        exp((mu[loc_idx[i]] - (0.5 * (sigma[loc_idx[i]]^2)) - (((log(lens[i]) - mu[loc_idx[i]] - log(rel_size[j]))^2) / (2 * sigma[loc_idx[i]]^2))))
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

  sig_o_mu <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_sigma <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_mu)
  ADREPORT(sig_o_sigma)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_mu", "eps_sigma"))
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

#### Bi-normal ####

pars <- list(
  ln_k1 = 5.25,
  ln_k2 = 2.5,
  ln_k3 = 5.5,
  ln_k4 = 3,
  ln_c = -.9,
  eps_k1 = rep(0, mesh$n),
  eps_k2 = rep(0, mesh$n),
  eps_k3 = rep(0, mesh$n),
  eps_k4 = rep(0, mesh$n),
  eps_c = rep(0, mesh$n),
  log_tau1 = 1,
  log_tau2 = 1,
  log_tau3 = 1,
  log_tau4 = 1,
  log_tau5 = 1,
  log_kappa1 = -1,
  log_kappa2 = -1,
  log_kappa3 = -1,
  log_kappa4 = -1,
  log_kappa5 = -1
)

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

  # range1 <- sqrt(8) / exp(log_kappa1)
  # range2 <- sqrt(8) / exp(log_kappa2)
  #
  # sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  # sig_o_k2 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))
  #
  # ADREPORT(range1)
  # ADREPORT(range2)
  # ADREPORT(sig_o_k1)
  # ADREPORT(sig_o_k2)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_k1", "eps_k2", "eps_k3", "eps_k4", "eps_c"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)

# Model does not converge


#### End of script ####
