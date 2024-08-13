# Code to produce mesh-specific selection curves and total selection curve
# for by-site binormal curve with multifilament data

library(RTMB)
library(tidyverse)

# Pull in data
by_net_multi <- read.csv("data/multi_catch_by_net.csv")
by_net_multi <- subset(by_net_multi, net_id != 2012005)
# remove this site since coordinates cannot be trusted

meshes <- seq(2, 5, 0.25) * 25.4 # stretch measure mesh size in mm
n_sites <- length(unique(by_net_multi$net_id)) # number of sites

# Pull in model estimates from RTMB
rep <- readRDS("data/rep_binorm_bysite_multi_all.rds")

# Lengths to plot

plot_lengths <- seq(100, 900, 0.1)

rel_size <- meshes / meshes[1] # mesh sizes relative to first panel

n_obs <- length(plot_lengths) * n_sites # number of observations = number of lengths to plot x number of sites

len_seq <- rep(plot_lengths, n_sites) # sequence of length classes repeated for each site


################################################
# add a mean selectivity curve using only global averages

ln_k1 <- unname(rep$par.fixed[1]) # global estimate for ln_k1. unname() removes TMB stuff

ln_k2 <- unname(rep$par.fixed[2]) # global estimate for ln_k2.

ln_k3 <- unname(rep$par.fixed[3])

ln_k4 <- unname(rep$par.fixed[4])

ln_c <- unname(rep$par.fixed[5])

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
  "panel_11", "panel_12", "panel_13"
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


plot_mean$r_tot <- rowSums(plot_mean[, 8:20]) # calculate total selectivity for mean curve

plot_mean$r_sel_max <- rep(max(plot_mean$r_tot), length(plot_lengths))

plot_mean$r_tot_scaled <- plot_mean$r_tot / max(plot_mean$r_tot) # scale total selectivity to one


plot_mean_multi <- plot_mean

################################################
plot_meshes_multi <- pivot_longer(plot_mean_multi,
  cols = panel_1:panel_13,
  names_to = "mesh",
  values_to = "r_sel"
)

plot_meshes_multi %>%
  group_by(mesh) %>%
  mutate(plot_sel = r_sel / max(r_sel)) %>%
  ungroup() -> plot_meshes_multi

plot_meshes_multi <- data.frame(plot_meshes_multi)

# Convert mesh to factor with desired levels
plot_meshes_multi$mesh <- factor(plot_meshes_multi$mesh,
  levels = c("panel_1", "panel_2", "panel_3", "panel_4", "panel_5", "panel_6", "panel_7", "panel_8", "panel_9", "panel_10", "panel_11", "panel_12", "panel_13"),
  labels = c("51", "57", "64", "70", "76", "83", "89", "95", "102", "108", "114", "121", "127")
)

p1 <- ggplot(
  plot_meshes_multi[plot_meshes_multi$site_id == "global", ],
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

p2 <- ggplot(data = plot_mean_multi, aes(x = length)) +
  geom_line(aes(y = r_tot_scaled), size = 1.5) +
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
