# Script to calculate and plot deviance residuals for top models for
# monofilament and multifilament catch data

#### Multifilament data ####

# Read in multifilament data
data <- read.csv("data/multi_catch_by_net.csv")

# need to remove 2012005 due to unreliable coordinates
data <- data[data$net_id != 2012005, ]

meshes <- seq(2.0, 5.0, by = 0.25) # range of stretch measure meshes in inches
meshes <- meshes * 25.4 # convert to mm

n_sites <- length(unique(data$net_id))
catches <- as.matrix(data[, which(grepl("mesh", colnames(data)))])
colnames(catches) <- NULL

# Create list for TMB
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

# Top model = non-spatial binormal function with all parameters varying

# Pull in model estimates from RTMB
rep <- readRDS("data/rep_binorm_bysite_multi_all.rds")
obj <- readRDS("data/obj_binorm_bysite_multi_all.rds")
opt <- readRDS("data/opt_binorm_bysite_multi_all.rds")

# Code modified from Millar's website
# https://www.stat.auckland.ac.nz/~millar/selectware/code.html

O <- data$catches

phi <- obj$report(obj$env$last.par.best)$phi_mat

E <- apply(O, 1, sum, na.rm = TRUE) * phi

wk <- O * log(O / E)
wk[is.na(wk)] <- 0

Dev.resids <- sign(O - E) * sqrt(2 * (E - O + wk))

plot_resids <- cbind(data$site, data$lens, Dev.resids)

colnames(plot_resids) <- c(
  "site_id",
  "length",
  "panel_1",
  "panel_2",
  "panel_3",
  "panel_4",
  "panel_5",
  "panel_6",
  "panel_7",
  "panel_8",
  "panel_9",
  "panel_10",
  "panel_11",
  "panel_12",
  "panel_13"
)

library(tidyverse)

plot_resids_long <- as.data.frame(plot_resids) %>%
  pivot_longer(!c("site_id", "length"),
    names_to = "panel", values_to = "dev_resid"
  )

plot_resids_long$sign <- ifelse(plot_resids_long$dev_resid >= 0, "positive", "negative")

plot_resids_long$site_id <- as.factor(plot_resids_long$site_id)

ggplot(
  plot_resids_long,
  aes(
    x = length,
    y = panel,
    size = dev_resid,
    shape = sign
  )
) +
  geom_point() +
  scale_size_area(max_size = 10) +
  scale_shape_manual(values = c(
    "positive" = 16, # 16 = black circle
    "negative" = 1
  )) + # 1 = open circle
  scale_y_discrete(
    name = "Panel size (stretch, mm)",
    limits = c(
      "panel_1", "panel_2", "panel_3", "panel_4",
      "panel_5", "panel_6", "panel_7", "panel_8",
      "panel_9", "panel_10", "panel_11", "panel_12",
      "panel_13"
    ), # reorder panels on axis
    labels = c(
      "panel_1" = "51", "panel_2" = "57", "panel_3" = "64",
      "panel_4" = "70", "panel_5" = "76", "panel_6" = "83",
      "panel_7" = "89", "panel_8" = "95", "panel_9" = "102",
      "panel_10" = "108", "panel_11" = "114", "panel_12" = "121",
      "panel_13" = "127"
    )
  ) + # relabel with sizes
  xlab("Length (mm)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black")
  ) +
  facet_wrap(. ~ site_id)


#### Monofilament data ####

library(RTMB)

# Read in monofilament data
data <- read.csv("data/by_net_mono.csv")

data <- data[data$CRNII != 2012005, ] # remove this site due to unreliable coordinates

meshes <- seq(1.5, 7.0, by = 0.5) # range of stretch measure meshes in inches
meshes <- meshes * 25.4 # convert to mm

n_sites <- length(unique(data$CRNII))
catches <- as.matrix(data[, -c(1:2, 15:16)])
colnames(catches) <- NULL

# Create list for TMB
data <- list(
  mesh = meshes,
  n_site = length(unique(data$CRNII)),
  n_obs = length(data[, 1]),
  site = rep(
    seq(
      from = 1,
      to = length(unique(data$CRNII)),
      length.out = length(unique(data$CRNII))
    ),
    each = length(unique(data$lbin))
  ),
  lens = data[, 2],
  catches = catches
)

colnames(data$catches) <- NULL
data$rel_size <- data$mesh / data$mesh[1]

# Pull in model estimates from RTMB
rep <- readRDS("data/rep_binorm_bysite_mono_all.rds")
obj <- readRDS("data/obj_binorm_bysite_mono_all.rds")
opt <- readRDS("data/opt_binorm_bysite_mono_all.rds")

# Code modified from Millar's website
# https://www.stat.auckland.ac.nz/~millar/selectware/code.html

O <- data$catches

phi <- obj$report(obj$env$last.par.best)$phi_mat

E <- apply(O, 1, sum, na.rm = TRUE) * phi

wk <- O * log(O / E)
wk[is.na(wk)] <- 0

Dev.resids <- sign(O - E) * sqrt(2 * (E - O + wk))

plot_resids <- cbind(data$site, data$lens, Dev.resids)

colnames(plot_resids) <- c(
  "site_id",
  "length",
  "panel_1",
  "panel_2",
  "panel_3",
  "panel_4",
  "panel_5",
  "panel_6",
  "panel_7",
  "panel_8",
  "panel_9",
  "panel_10",
  "panel_11",
  "panel_12"
)

library(tidyverse)

plot_resids_long <- as.data.frame(plot_resids) %>%
  pivot_longer(!c("site_id", "length"),
    names_to = "panel", values_to = "dev_resid"
  )

plot_resids_long$sign <- ifelse(plot_resids_long$dev_resid >= 0, "positive", "negative")

plot_resids_long$site_id <- as.factor(plot_resids_long$site_id)

ggplot(
  plot_resids_long,
  aes(
    x = length,
    y = panel,
    size = dev_resid,
    shape = sign
  )
) +
  geom_point() +
  scale_size_area(max_size = 10) +
  scale_shape_manual(values = c(
    "positive" = 16, # 16 = black circle
    "negative" = 1
  )) + # 1 = open circle
  scale_y_discrete(
    name = "Panel size (stretch, mm)",
    limits = c(
      "panel_1", "panel_2", "panel_3", "panel_4",
      "panel_5", "panel_6", "panel_7", "panel_8",
      "panel_9", "panel_10", "panel_11", "panel_12"
    ), # reorder panels on axis
    labels = c(
      "panel_1" = "38", "panel_2" = "51", "panel_3" = "64",
      "panel_4" = "76", "panel_5" = "89", "panel_6" = "102",
      "panel_7" = "114", "panel_8" = "127", "panel_9" = "140",
      "panel_10" = "152", "panel_11" = "165", "panel_12" = "178"
    )
  ) + # relabel with sizes
  xlab("Length (mm)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black")
  ) +
  facet_wrap(. ~ site_id)


#### End of script ####
