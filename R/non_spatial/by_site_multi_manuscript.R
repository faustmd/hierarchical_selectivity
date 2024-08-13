# Fit hierarchical selectivity functions to by-site catch data from
# multifilament gillnet catch data from gear comparison study 2007-2008 & 2011-2013.
# Script currently allows all parameters to vary.

library(RTMB)

# Read in multifilament data
data <- read.csv("data/multi_catch_by_net.csv")

# need to remove 2012005 due to unreliable coordinates
data <- data[data$net_id != 2012005, ]

meshes <- seq(2.0, 5.0, by = 0.25) # range of stretch measure meshes in inches
meshes <- meshes * 25.4 # convert to mm

n_sites <- length(unique(data$net_id))
catches <- as.matrix(data[, which(grepl("mesh", colnames(data)))])
colnames(catches) <- NULL

# Create list for RTMB
data <- list(
  mesh = meshes,
  n_site = length(unique(data$net_id)),
  n_obs = length(data[, 1]),
  site = rep(
    seq(
      from = 1,
      to = length(unique(data$net_id)),
      length.out = length(unique(data$net_id))
    ),
    each = length(unique(data$length_class))
  ),
  lens = data[, 2],
  catches = catches
)

colnames(data$catches) <- NULL
data$rel_size <- data$mesh / data$mesh[1]

### Normal location ####

# Normal curve, but assumes one variance among all meshes rather than variance
# scaling with mesh size

pars <- list(
  ln_k = log(225),
  ln_sigma = log(75),
  k_dev = rep(0, data$n_site),
  sigma_dev = rep(0, data$n_site),
  ln_sd_k = log(1),
  ln_sd_sigma = log(1)
)

f <- function(pars) {
  getAll(data, pars)
  catches <- OBS(catches)
  jnll <- 0
  jnll <- jnll - sum(dnorm(k_dev, 0, exp(ln_sd_k), TRUE))
  jnll <- jnll - sum(dnorm(sigma_dev, 0, exp(ln_sd_sigma), TRUE))
  k <- exp(ln_k + k_dev)
  sigma <- exp(ln_sigma + sigma_dev)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k[site[i]] * rel_size[j])^2 /
        (2 * sigma[site[i]]^2))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  ADREPORT(phi_mat)
  jnll
}

# Fit the model
obj <- MakeADFun(f, pars, random = c("k_dev", "sigma_dev"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$objective

# run with Newton steps to aid convergence
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

norm_pars <- list(
  ln_k1 = log(225),
  ln_k2 = log(75),
  k1_dev = rep(0, data$n_site),
  k2_dev = rep(0, data$n_site),
  ln_sd_k1 = log(1),
  ln_sd_k2 = log(1)
)

f_norm <- function(pars) {
  getAll(data, pars)
  catches <- OBS(catches)
  jnll <- 0
  jnll <- jnll - sum(dnorm(k1_dev, 0, exp(ln_sd_k1), TRUE))
  jnll <- jnll - sum(dnorm(k2_dev, 0, exp(ln_sd_k2), TRUE))
  k1 <- exp(ln_k1 + k1_dev)
  k2 <- exp(ln_k2 + k2_dev)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[site[i]] * rel_size[j])^2 /
        (2 * k2[site[i]]^2 * rel_size[j]^2))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(phi_mat)
  jnll
}

# Fit model
obj <- MakeADFun(f_norm, norm_pars, random = c("k1_dev", "k2_dev"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

# run with Newton steps to aid convergence
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

gamma_pars <- list(
  ln_alpha = 4,
  ln_k = 2,
  alpha_dev = rep(0, n_sites),
  k_dev = rep(0, n_sites),
  ln_sd_alpha = 1,
  ln_sd_k = 1
)

f_gamma <- function(pars) {
  getAll(data, pars)
  catches <- OBS(catches)
  jnll <- 0
  jnll <- jnll - sum(dnorm(alpha_dev, 0, exp(ln_sd_alpha), TRUE))
  jnll <- jnll - sum(dnorm(k_dev, 0, exp(ln_sd_k), TRUE))
  alpha <- exp(ln_alpha + alpha_dev)
  k <- exp(ln_k + k_dev)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- ((lens[i] / ((alpha[site[i]] - 1) * k[site[i]] * rel_size[j]))^(alpha[site[i]] - 1)) *
        (exp(alpha[site[i]] - 1 - (lens[i] / (k[site[i]] * rel_size[j]))))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(phi_mat)
  jnll
}

# Fit the model
obj <- MakeADFun(f_gamma, gamma_pars, random = c("alpha_dev", "k_dev"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$objective

# run with Newton steps to aid convergence
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

#### Log-normal ####

lognorm_pars <- list(
  ln_mu = 5.6,
  ln_sigma = 0.2,
  mu_dev = rep(0, n_sites),
  sigma_dev = rep(0, n_sites),
  ln_sd_mu = -1,
  ln_sd_sigma = -1
)

f_lognorm <- function(pars) {
  getAll(data, pars)
  catches <- OBS(catches)
  jnll <- 0
  jnll <- jnll - sum(dnorm(mu_dev, 0, exp(ln_sd_mu), TRUE))
  jnll <- jnll - sum(dnorm(sigma_dev, 0, exp(ln_sd_sigma), TRUE))
  mu <- (ln_mu + mu_dev)
  sigma <- (ln_sigma + sigma_dev)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- (rel_size[j] / lens[i]) *
        exp((mu[site[i]] - (0.5 * (sigma[site[i]]^2)) - (((log(lens[i]) - mu[site[i]] - log(rel_size[j]))^2) / (2 * sigma[site[i]]^2))))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(phi_mat)
  jnll
}

obj <- MakeADFun(f_lognorm, lognorm_pars, random = c("mu_dev", "sigma_dev"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

# run with Newton steps to aid convergence
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

#### Bi-normal ####

binorm_pars <- list(
  ln_k1 = 5.5,
  ln_k2 = 3.2,
  ln_k3 = 5.8,
  ln_k4 = 4.3,
  ln_c = -0.3,
  k1_dev = rep(0, n_sites),
  k2_dev = rep(0, n_sites),
  k3_dev = rep(0, n_sites),
  k4_dev = rep(0, n_sites),
  c_dev = rep(0, n_sites),
  ln_sd_k1 = 0.1,
  ln_sd_k2 = 0.1,
  ln_sd_k3 = 0.1,
  ln_sd_k4 = 0.1,
  ln_sd_c = 0.1
)

f_binorm <- function(pars) {
  getAll(data, pars)
  catches <- OBS(catches)
  jnll <- 0
  jnll <- jnll - sum(dnorm(k1_dev, 0, exp(ln_sd_k1), TRUE))
  jnll <- jnll - sum(dnorm(k2_dev, 0, exp(ln_sd_k2), TRUE))
  jnll <- jnll - sum(dnorm(k3_dev, 0, exp(ln_sd_k3), TRUE))
  jnll <- jnll - sum(dnorm(k4_dev, 0, exp(ln_sd_k4), TRUE))
  jnll <- jnll - sum(dnorm(c_dev, 0, exp(ln_sd_c), TRUE))
  k1 <- exp(ln_k1 + k1_dev)
  k2 <- exp(ln_k2 + k2_dev)
  k3 <- exp(ln_k3 + k3_dev)
  k4 <- exp(ln_k4 + k4_dev)
  c <- exp(ln_c + c_dev)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[site[i]] * rel_size[j])^2 /
        (2 * k2[site[i]]^2 * rel_size[j]^2)) +
        (c[site[i]] * exp(-(lens[i] - k3[site[i]] * rel_size[j])^2 /
          (2 * k4[site[i]]^2 * rel_size[j]^2)))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(phi_mat)
  jnll
}

# Fit the model
obj <- MakeADFun(f_binorm, binorm_pars, random = c("k1_dev", "k2_dev", "k3_dev", "k4_dev", "c_dev"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

# run with Newton steps to aid convergence
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

#### End of script ####
