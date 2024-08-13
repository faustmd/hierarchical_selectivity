# Fit selectivity functions to pooled catch data for multifilament
# catch data from gear comparison study 2007-2008 & 2011-2013.
# Script represents approach outlined by Millar and Fryer (1999)
# to estimating gillnet selectivity.
#
# These parameters are used to calculate log-likelihood for
# set-specific data, which allows apples to apples comparison via
# model selection with non-spatial & spatial versions of these functions.

library(RTMB)

# Pull in pooled catch data for monofilament gear
data_short <- read.csv("data/mono_catch_pooled.csv")

meshes <- seq(1.5, 7.0, by = 0.5) # range of stretch measure meshes in inches
meshes <- meshes * 25.4 # convert to mm

# Create list for RTMB
data_pooled <- list(
  mesh = meshes,
  lens = data_short$lbin,
  catches = as.matrix(data_short[, -c(1)])
)
colnames(data_pooled$catches) <- NULL
data_pooled$rel_size <- data_pooled$mesh / data_pooled$mesh[1]

# By-site monofilament catch data to use with pooled parameter estimates
# this is used in model selection process
data_long <- read.csv("data/mono_catch_by_net.csv")

# remove this net due to wrong coordinates
data_long <- data_long[data_long$CRNII != 2012005, ]
data_long <- data_long[, -c(15, 16)] # removes lat and lon columns

# Create list for TMB
data_by_net <- list(
  mesh = meshes,
  lens = data_long$lbin,
  catches = as.matrix(data_long[, -c(1:2)])
)
colnames(data_by_net$catches) <- NULL
data_by_net$rel_size <- data_by_net$mesh / data_by_net$mesh[1]

#### Pooled normal location ####

# This block of code fits model to pooled (in Millar's way) catch data

pars <- list(
  ln_k = log(225),
  ln_sigma = log(75)
)

f <- function(pars) {
  getAll(data_pooled, pars)
  catches <- OBS(catches)
  jnll <- 0
  k <- exp(ln_k)
  sigma <- exp(ln_sigma)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k * rel_size[j])^2 /
        (2 * sigma^2))
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
  REPORT(jnll) # to calculate log-likelihood with Millar pooled parameters
  jnll
}

# Fit the model using TMB
obj <- MakeADFun(f, pars)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
rep <- sdreport(obj)
rep

# Now use these parameter estimates to calculate log-likelihood for by-net
# version of the monofilament catch data

pars <- list(
  ln_k = rep$par.fixed[1],
  ln_sigma = rep$par.fixed[2]
)

f <- function(pars) {
  getAll(data_by_net, pars)
  catches <- OBS(catches)
  jnll <- 0
  k <- exp(ln_k)
  sigma <- exp(ln_sigma)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k * rel_size[j])^2 /
        (2 * sigma^2))
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
  REPORT(jnll) # to calculate log-likelihood with Millar pooled parameters
  jnll
}

# calculate log likelihood
obj <- MakeADFun(f, pars)
obj$report()$jnll

#### Pooled normal ####

# First fit model following Millar by pooling all catch data

# starting values for parameters
pars <- list(
  ln_k1 = log(250),
  ln_k2 = log(50)
)

f <- function(pars) {
  getAll(data_pooled, pars)
  catches <- OBS(catches)
  jnll <- 0
  k1 <- exp(ln_k1)
  k2 <- exp(ln_k2)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1 * rel_size[j])^2 /
        (2 * k2^2 * rel_size[j]^2))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(jnll)
  jnll
}

# Fit the model using TMB
obj <- MakeADFun(f, pars)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
rep <- sdreport(obj)
rep

# Now use these parameter estimates to calculate log-likelihood for by-net
# version of the monofilament catch data

# starting values for parameters
pars <- list(
  ln_k1 = rep$par.fixed[1],
  ln_k2 = rep$par.fixed[2]
)

f <- function(pars) {
  getAll(data_by_net, pars) # data_by_net is only change
  catches <- OBS(catches)
  jnll <- 0
  k1 <- exp(ln_k1)
  k2 <- exp(ln_k2)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1 * rel_size[j])^2 /
        (2 * k2^2 * rel_size[j]^2))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(jnll)
  jnll
}

# Now calculate log-likelihood for by-net data given pooled estimates above
obj <- MakeADFun(f, pars)
obj$report()$jnll

#### Pooled log-normal ####

# Starting parameters
pars <- list(
  ln_mu = 5.5,
  ln_sigma = 1
)

f_lognorm <- function(pars) {
  getAll(data_pooled, pars) # data_pooled is only change
  catches <- OBS(catches)
  jnll <- 0
  mu <- (ln_mu)
  sigma <- (ln_sigma)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- (rel_size[j] / lens[i]) *
        exp((mu - (0.5 * (sigma^2)) - (((log(lens[i]) - mu - log(rel_size[j]))^2) / (2 * sigma^2))))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(jnll)
  jnll
}

# Fit model using TMB
obj <- MakeADFun(f_lognorm, pars)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
rep <- sdreport(obj)
rep

# Now use these parameter estimates to calculate log-likelihood for by-net
# version of the monofilament catch data

# Starting parameters
pars <- list(
  ln_mu = rep$par.fixed[1],
  ln_sigma = rep$par.fixed[2]
)

f_lognorm <- function(pars) {
  getAll(data_by_net, pars) # data_by_net is only change
  catches <- OBS(catches)
  jnll <- 0
  mu <- (ln_mu)
  sigma <- (ln_sigma)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- (rel_size[j] / lens[i]) *
        exp((mu - (0.5 * (sigma^2)) - (((log(lens[i]) - mu - log(rel_size[j]))^2) / (2 * sigma^2))))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(jnll)
  jnll
}

# Calculate log-likelihood using pooled parameters above, but with by-net data
obj <- MakeADFun(f_lognorm, pars)
obj$report()$jnll

#### Pooled gamma ####

# First fit model to pooled data following Millar examples

# Starting parameters
pars <- list(
  ln_alpha = 4,
  ln_k = 2
)

f_gamma <- function(pars) {
  getAll(data_pooled, pars)
  catches <- OBS(catches)
  jnll <- 0
  alpha <- exp(ln_alpha)
  k <- exp(ln_k)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- ((lens[i] / ((alpha - 1) * k * rel_size[j]))^(alpha - 1)) *
        (exp(alpha - 1 - (lens[i] / (k * rel_size[j]))))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(jnll)
  jnll
}

obj <- MakeADFun(f_gamma, pars)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
rep <- sdreport(obj)
rep

# Now calculate log-likelihood using these (pooled) parameters with by-net data

# Starting parameters
pars <- list(
  ln_alpha = rep$par.fixed[1],
  ln_k = rep$par.fixed[2]
)

f_gamma <- function(pars) {
  getAll(data_by_net, pars) # Only change is data
  catches <- OBS(catches)
  jnll <- 0
  alpha <- exp(ln_alpha)
  k <- exp(ln_k)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- ((lens[i] / ((alpha - 1) * k * rel_size[j]))^(alpha - 1)) *
        (exp(alpha - 1 - (lens[i] / (k * rel_size[j]))))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE))
  REPORT(jnll)
  jnll
}

# calculate log-likelihood using pooled parameters
obj <- MakeADFun(f_gamma, pars)
obj$report()$jnll


#### Pooled bi-normal ####

# Fit model to pooled data following Millar examples

pars <- list(
  ln_k1 = 5.5,
  ln_k2 = 2.5,
  ln_k3 = 5.9,
  ln_k4 = 3.5,
  ln_c = 0.5
)

f_binorm <- function(pars) {
  getAll(data_pooled, pars)
  catches <- OBS(catches)
  jnll <- 0
  k1 <- exp(ln_k1)
  k2 <- exp(ln_k2)
  k3 <- exp(ln_k3)
  k4 <- exp(ln_k4)
  c <- exp(ln_c)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1 * rel_size[j])^2 /
        (2 * k2^2 * rel_size[j]^2)) +
        (c * exp(-(lens[i] - k3 * rel_size[j])^2 /
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
  REPORT(jnll)
  jnll
}

obj <- MakeADFun(f_binorm, pars)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
rep <- sdreport(obj)
rep

# Now calculate log-likelihood for by-net data using pooled parameters

pars <- list(
  ln_k1 = rep$par.fixed[1],
  ln_k2 = rep$par.fixed[2],
  ln_k3 = rep$par.fixed[3],
  ln_k4 = rep$par.fixed[4],
  ln_c = rep$par.fixed[5]
)

f_binorm <- function(pars) {
  getAll(data_by_net, pars) # only change is to what data are used
  catches <- OBS(catches)
  jnll <- 0
  k1 <- exp(ln_k1)
  k2 <- exp(ln_k2)
  k3 <- exp(ln_k3)
  k4 <- exp(ln_k4)
  c <- exp(ln_c)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1 * rel_size[j])^2 /
        (2 * k2^2 * rel_size[j]^2)) +
        (c * exp(-(lens[i] - k3 * rel_size[j])^2 /
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
  REPORT(jnll)
  jnll
}

# Calculate log-likelihood
obj <- MakeADFun(f_binorm, pars)
obj$report()$jnll

#### End of script ####
