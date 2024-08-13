# Create a plot showing mesh-specific selection curves and total selection curve
# across all panels for top model (non-spatial bi-normal) for monofilament net

library(RTMB)
library(tidyverse)

by_net_mono <- read.csv("data/by_net_mono.csv")

by_net_mono <- by_net_mono[by_net_mono$CRNII != 2012005, ]

meshes <- seq(1.5, 7.0, by = 0.5)
meshes <- meshes * 25.4

n_sites <- length(unique(by_net_mono$CRNII)) # number of sites

# Bi-normal curve

# pull in RTMB estimates

rep <- readRDS("data/rep_binorm_bysite_mono_all.rds")

# lengths to plot

plot_lengths <- seq(100, 900, 0.1)

rel_size <- meshes / meshes[1] # mesh sizes relative to first panel

n_obs <- length(plot_lengths) * n_sites # number of observations = number of lengths to plot x number of sites

len_seq <- rep(plot_lengths, n_sites) # sequence of length classes repeated for each site

site_id <- rep(unique(by_net_mono$CRNII), each = length(plot_lengths))

# parameter estimates from by-site bi-normal curve

ln_k1 <- unname(rep$par.fixed[1]) # global estimate for ln_k1. unname() removes TMB stuff

ln_k2 <- unname(rep$par.fixed[2]) # global estimate for ln_k2.

ln_k3 <- unname(rep$par.fixed[3])

ln_k4 <- unname(rep$par.fixed[4])

ln_c <- unname(rep$par.fixed[5])

k1_devs <- unname(rep$par.random[names(rep$par.random) == "k1_dev"]) # site-specific deviations for sites 1-63

k2_devs <- unname(rep$par.random[names(rep$par.random) == "k2_dev"]) # site-specific deviations for sites 1:63

k3_devs <- unname(rep$par.random[names(rep$par.random) == "k3_dev"])

k4_devs <- unname(rep$par.random[names(rep$par.random) == "k4_dev"])

c_devs <- unname(rep$par.random[names(rep$par.random) == "c_dev"])

k1_net <- exp(ln_k1 + k1_devs) # back-transformed alpha for each site

k1_seq <- rep(k1_net, each = length(plot_lengths)) # sequence of site-specific k1 for plotting

k2_net <- exp(ln_k2 + k2_devs) # site-specific k2 estimates

k2_seq <- rep(k2_net, each = length(plot_lengths)) # sequence of site-specific k2 for plotting

k3_net <- exp(ln_k3 + k3_devs)

k3_seq <- rep(k3_net, each = length(plot_lengths))

k4_net <- exp(ln_k4 + k4_devs)

k4_seq <- rep(k4_net, each = length(plot_lengths))

c_net <- exp(ln_c + c_devs)

c_seq <- rep(c_net, each = length(plot_lengths))

mesh_mat <- matrix(data = NA, nrow = n_obs, ncol = length(meshes)) # matrix of length n_obs and width = number of meshed - just to add to data frame below

plot_data <- data.frame(cbind(
  len_seq, site_id, k1_seq,
  k2_seq, k3_seq, k4_seq,
  c_seq, mesh_mat
)) # pull elements into a single object

colnames(plot_data) <- c(
  "length", "site_id", "k1_net", "k2_net", "k3_net", "k4_net",
  "c_net", "panel_1", "panel_2", "panel_3", "panel_4", "panel_5",
  "panel_6", "panel_7", "panel_8", "panel_9", "panel_10", "panel_11",
  "panel_12"
)

for (i in 1:length(rel_size)) {
  plot_data[, i + 7] <- exp(-((plot_data$length - (plot_data$k1_net * rel_size[i]))^2) / (2 * (plot_data$k2_net^2) * (rel_size[i]^2))) + plot_data$c_net * exp(-((plot_data$length - (plot_data$k3_net * rel_size[i]))^2) / (2 * (plot_data$k4_net^2) * (rel_size[i]^2)))
}

plot_data$r_tot <- rowSums(plot_data[, 8:19]) # calculate total selectivity for a given size by each net

max_sel <- plot_data %>%
  group_by(site_id) %>%
  summarize(max_sel = max(r_tot)) # pull out maximum selectivity by site

max_sel <- unname(max_sel$max_sel) # create vector of max selectivities by site

plot_data$r_sel_max <- rep(max_sel, each = length(plot_lengths)) # create sequence by site repeating over lengths

plot_data$r_tot_scaled <- plot_data$r_tot / plot_data$r_sel_max # scale total selectivity to one for each site

# add a mean selectivity curve using only global averages

k1_mean <- exp(ln_k1)

k2_mean <- exp(ln_k2)

k3_mean <- exp(ln_k3)

k4_mean <- exp(ln_k4)

c_mean <- exp(ln_c)

mesh_mat <- matrix(data = NA, nrow = length(plot_lengths), ncol = length(rel_size))

plot_mean <- data.frame(cbind(
  plot_lengths,
  rep("global", length(plot_lengths)),
  rep(k1_mean, length(plot_lengths)),
  rep(k2_mean, length(plot_lengths)),
  rep(k3_mean, length(plot_lengths)),
  rep(k4_mean, length(plot_lengths)),
  rep(c_mean, length(plot_lengths)),
  mesh_mat
))

colnames(plot_mean) <- c(
  "length", "site_id", "k1_net", "k2_net", "k3_net", "k4_net",
  "c_net", "panel_1", "panel_2", "panel_3", "panel_4", "panel_5",
  "panel_6", "panel_7", "panel_8", "panel_9", "panel_10",
  "panel_11", "panel_12"
)

# below are characters so need to convert to numeric

plot_mean$length <- as.numeric(plot_mean$length)
plot_mean$k1_net <- as.numeric(plot_mean$k1_net)
plot_mean$k2_net <- as.numeric(plot_mean$k2_net)
plot_mean$k3_net <- as.numeric(plot_mean$k3_net)
plot_mean$k4_net <- as.numeric(plot_mean$k4_net)
plot_mean$c_net <- as.numeric(plot_mean$c_net)

# calculate mean selectivity
for (i in 1:length(rel_size)) {
  plot_mean[, i + 7] <- exp(-((plot_mean$length - (plot_mean$k1_net * rel_size[i]))^2) / (2 * (plot_mean$k2_net^2) * (rel_size[i]^2))) + plot_mean$c_net * exp(-((plot_mean$length - (plot_mean$k3_net * rel_size[i]))^2) / (2 * (plot_mean$k4_net^2) * (rel_size[i]^2)))
}


plot_mean$r_tot <- rowSums(plot_mean[, 8:19]) # calculate total selectivity for mean curve

plot_mean$r_sel_max <- rep(max(plot_mean$r_tot), length(plot_lengths))

plot_mean$r_tot_scaled <- plot_mean$r_tot / max(plot_mean$r_tot) # scale total selectivity to one

# now combine plot_mean with plot_data to have it all in one spot

plot_data <- rbind(plot_mean, plot_data) # gotta fix names so they match - I outsmarted myself - chance alpha_mean to alpha_net, etc.

plot_data <- select(plot_data, 1:7, 22)

plot_mean_mono <- plot_mean

plot_data_mono <- plot_data

plot_meshes_mono <- pivot_longer(plot_mean_mono,
  cols = panel_1:panel_12,
  names_to = "mesh",
  values_to = "r_sel"
)
plot_meshes_mono %>%
  group_by(mesh) %>%
  mutate(plot_sel = r_sel / max(r_sel)) %>%
  ungroup() -> plot_meshes_mono

plot_meshes_mono <- data.frame(plot_meshes_mono)

# Convert mesh to factor with desired levels
plot_meshes_mono$mesh <- factor(plot_meshes_mono$mesh,
  levels = c("panel_1", "panel_2", "panel_3", "panel_4", "panel_5", "panel_6", "panel_7", "panel_8", "panel_9", "panel_10", "panel_11", "panel_12"),
  labels = c("38", "51", "64", "76", "89", "102", "114", "127", "140", "152", "165", "179")
)

p1 <- ggplot(
  plot_meshes_mono[plot_meshes_mono$site_id == "global", ],
  aes(x = length, y = plot_sel, color = mesh)
) +
  geom_line(size = 1.5) +
  theme_classic() +
  labs(
    x = "Total length (mm)", y = "Relative selectivity",
    color = expression(atop(
      "Mesh size",
      paste("(stretch measure; mm)")
    ))
  ) +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1),
    legend.title = element_text(color = "black", size = 18),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top"
  )
p1
p2 <- ggplot(data = plot_data_mono, aes(x = length)) +
  geom_line(
    data = plot_data_mono[plot_data_mono$site_id == "global", ],
    aes(y = r_tot_scaled), size = 1.5
  ) +
  theme_classic() +
  labs(x = "Total length (mm)", y = "Relative selectivity") +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1)
  )
p2
library(patchwork)

p1 / p2

#### End of script ####
