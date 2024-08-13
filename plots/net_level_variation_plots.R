# Code to produce figures showing net-specific information gained from
# hierarchical approach to indirect estimation of contact selectivity for
# monofilament and multifilament gill nets.
# Code produces figures 7, 8, and 9 in submitted manuscript.

library(RTMB)
library(tidyverse)
library(dplyr)

# First figure shows length at peak selectivity (relative selectivity = 1.0)
# for each net and gear type

# Bring in .rds objects for multi and mono data with predicted selectivities
# at lengths for each net

# Multifilament data
sel_multi <- readRDS(file = "data/by_site_selectivity_multi.RDS")

# Subset down to length at peak selectivity (i.e., r_tot_scaled = 1.0)
peak_len_multi <- subset(sel_multi, r_tot_scaled == max(r_tot_scaled))
peak_len_multi$gear <- rep("multi", length(peak_len_multi$site_id))

# Monofilament data
sel_mono <- readRDS(file = "data/by_site_selectivity_mono.RDS")

# Subset down to length at peak selectivity (i.e., r_tot_scaled = 1.0)
peak_len_mono <- subset(sel_mono, r_tot_scaled == max(r_tot_scaled))
peak_len_mono$gear <- rep("mono", length(peak_len_mono$site_id))

# Merge multi and mono data together for plotting
peak_len <- rbind(peak_len_multi[, c(1:2, 24)], peak_len_mono[, c(1:2, 23)])

#### First plot ####

# Make plot - should this be a violin plot like the next figure?

ggplot(data = peak_len) +
  theme_classic() +
  geom_violin(aes(x = gear, y = length), linewidth = 1.5) +
  # geom_boxplot(aes(x = gear, y = length), linewidth = 1.25)+
  geom_jitter(aes(x = gear, y = length, alpha = 0.5)) +
  xlab(NULL) +
  ylab("Total length (mm)") +
  scale_x_discrete(
    name = NULL,
    labels = c(
      "mono" = "Monofilament",
      "multi" = "Multifilament"
    )
  ) +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )

#### Second figure ####

# Shows distribution of selectivities at select lengths for both gear types

# Subset each gear's data down to selectivity at 300, 400, 500, and 600 mm
plot_multi <- sel_multi[sel_multi$length == 300 |
  sel_multi$length == 400 |
  sel_multi$length == 500 |
  sel_multi$length == 600 |
  sel_multi$length == 700, ]
plot_multi$gear <- "multi" # add gear label

plot_mono <- sel_mono[sel_mono$length == 300 |
  sel_mono$length == 400 |
  sel_mono$length == 500 |
  sel_mono$length == 600 |
  sel_mono$length == 700, ]
plot_mono$gear <- "mono" # add gear label

# combine data together for plotting
plot_gears <- rbind(plot_multi[, c(1:2, 23:24)], plot_mono[, c(1:2, 22:23)])

# Make the plot

ggplot(plot_gears) +
  geom_boxplot(aes(x = as.factor(length), y = r_tot_scaled, color = gear),
    outlier.shape = NA, linewidth = 1.25
  ) +
  geom_jitter(aes(x = as.factor(length), y = r_tot_scaled, color = gear),
    position = position_jitterdodge(), alpha = 0.5
  ) +
  theme_classic() +
  xlab("Length (mm)") +
  ylab("Relative selectivity") +
  scale_color_discrete(
    name = NULL,
    labels = c(
      "multi" = "Multifilament",
      "mono" = "Monofilament"
    )
  ) +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top",
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )

# Alternative with violin plot

ggplot(plot_gears) +
  geom_violin(aes(x = as.factor(length), y = r_tot_scaled, fill = gear),
    alpha = 0.5, linewidth = 1.25
  ) +
  geom_jitter(aes(x = as.factor(length), y = r_tot_scaled, color = gear),
    position = position_jitterdodge(), alpha = 0.5
  ) +
  theme_classic() +
  xlab("Total length (mm)") +
  ylab("Relative selectivity") +
  scale_color_discrete(
    name = NULL,
    labels = c(
      "mono" = "Monofilament",
      "multi" = "Multifilament"
    )
  ) +
  scale_fill_discrete(
    name = NULL,
    labels = c(
      "mono" = "Monofilament",
      "multi" = "Multifilament"
    )
  ) +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top",
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )

#### Third plot ####

# Plots to illustrate correlations among parameter deviations and net-level
# abiotic and biotic data

# Need to standardize covariates to aid in interpretation of plots!!!

# Pull in net-level covariates for plotting
covariates <- read.csv("data/plot_covariates.csv")
str(covariates)

covariates$SECCHI <- as.numeric(covariates$SECCHI)

# Center and scale covariates to aid interpreation on plots
covariates$tot_catch_s <- as.numeric(scale(covariates$tot_catch, center = TRUE, scale = TRUE))
covariates$wae_catch_s <- as.numeric(scale(covariates$wae_catch, center = TRUE, scale = TRUE))
covariates$depth_s <- as.numeric(scale(covariates$DEPTH, center = TRUE, scale = TRUE))
covariates$secchi_s <- as.numeric(scale(covariates$SECCHI, center = TRUE, scale = TRUE))

# Model output for bi-normal selectivity function for multifilament data
rep <- readRDS("data/rep_binorm_bysite_multi_all.rds")
obj <- readRDS("data/obj_binorm_bysite_multi_all.rds")
opt <- readRDS("data/opt_binorm_bysite_multi_all.rds")

plot_multi <- covariates[covariates$GEAR == 7, ]

plot_multi <- plot_multi[plot_multi$CRNII != 2012005 & plot_multi$CRNII != 2012018, ] # remove this curve since coordinates unreliable

# Parameter deviations for each parameter in selectivity function
plot_multi$k1_dev <- rep$par.random[names(rep$par.random) == "k1_dev"]
plot_multi$k2_dev <- rep$par.random[names(rep$par.random) == "k2_dev"]
plot_multi$k3_dev <- rep$par.random[names(rep$par.random) == "k3_dev"]
plot_multi$k4_dev <- rep$par.random[names(rep$par.random) == "k4_dev"]
plot_multi$c_dev <- rep$par.random[names(rep$par.random) == "c_dev"]

# This is only using standardized covariates
plot_multi_covariates <- plot_multi %>%
  pivot_longer(
    cols = c(tot_catch_s:secchi_s),
    names_to = "covariate",
    values_to = "covariate_value"
  )

plot_multi_pars <- plot_multi %>%
  pivot_longer(
    cols = c(k1_dev, k2_dev, k3_dev, k4_dev, c_dev),
    names_to = "parameter",
    values_to = "parameter_dev"
  )

# Merge the two pivoted data frames
plot_multi_long <- merge(plot_multi_covariates, plot_multi_pars, by = "CRNII")

# Select only columns worth keeping
plot_multi_long <- select(plot_multi_long, CRNII, covariate:covariate_value, parameter:parameter_dev)

plot_multi_long$gear <- "multi"

# Calculate Pearson correlation between standardized covariates
# and parameter deviations.

plot_multi <- plot_multi[complete.cases(plot_multi), ]

# Select columns ending with "_dev" and "_s"
dev_columns <- names(plot_multi)[grep("_dev$", names(plot_multi))]
s_columns <- names(plot_multi)[grep("_s$", names(plot_multi))]

# Initialize an empty data frame to store results
correlation_results <- data.frame()

# Loop through each combination of dev and s columns
for (dev_col in dev_columns) {
  for (s_col in s_columns) {
    # Calculate Pearson correlation
    correlation <- cor(plot_multi[[dev_col]], plot_multi[[s_col]])

    # Calculate Fisher transformation
    fisher_transform <- 0.5 * log((1 + correlation) / (1 - correlation))

    # Compute standard error
    n <- nrow(plot_multi)
    standard_error <- 1 / sqrt(n - 3)

    # Compute z value for 95% confidence interval
    z_value <- qnorm(0.975)

    # Compute margin of error
    margin_of_error <- z_value * standard_error

    # Compute lower and upper bounds of the confidence interval for the Fisher transformation
    lower_bound <- fisher_transform - margin_of_error
    upper_bound <- fisher_transform + margin_of_error

    # Apply inverse Fisher transformation to get the confidence interval for the correlation coefficient
    lower_ci <- (exp(2 * lower_bound) - 1) / (exp(2 * lower_bound) + 1)
    upper_ci <- (exp(2 * upper_bound) - 1) / (exp(2 * upper_bound) + 1)

    # Extract gear and parameter names
    gear <- plot_multi$GEAR[1] # Assuming all rows have the same gear
    parameter <- sub("_dev", "", dev_col)
    covariate <- sub("_s", "", s_col)

    # Create a data frame for the current correlation with confidence intervals
    correlation_df <- data.frame(
      gear = gear,
      parameter = parameter,
      covariate = covariate,
      correlation = correlation,
      lower_ci = lower_ci,
      upper_ci = upper_ci
    )

    # Append to the main results data frame
    correlation_results <- rbind(correlation_results, correlation_df)
  }
}

# Print the results
print(correlation_results)

# For facet labels
facet_labels <- list(
  covariate = c(
    "DEPTH" = "Lake depth (m)",
    "tot_catch" = "Total fishes caught",
    "wae_catch" = "Total walleye caught",
    "SECCHI" = "Secchi depth (m)"
  ),
  parameter = c(
    "k1_dev" = "k1",
    "k2_dev" = "k2",
    "k3_dev" = "k3",
    "k4_dev" = "k4",
    "c_dev" = "c"
  )
)

# If using standardized covariates
facet_labels <- list(
  covariate = c(
    "depth_s" = "Lake depth (m)",
    "tot_catch_s" = "Total fishes caught",
    "wae_catch_s" = "Total walleye caught",
    "secchi_s" = "Secchi depth (m)"
  ),
  parameter = c(
    "k1_dev" = "k1",
    "k2_dev" = "k2",
    "k3_dev" = "k3",
    "k4_dev" = "k4",
    "c_dev" = "c"
  )
)

# Filter out rows with NA values
plot_multi_filtered <- plot_multi_long[complete.cases(plot_multi_long), ]

# Standardized covariates and correlation labels

# Merge correlation estimates and uncertainty around it to the plotting data

plot_multi_filtered$parameter <- sub("_dev.*", "", plot_multi_filtered$parameter)
plot_multi_filtered$covariate <- sub("_s.*", "", plot_multi_filtered$covariate)

correlation_results$gear <- "multi"

merged_data <- merge(plot_multi_filtered,
  correlation_results,
  by = c("gear", "parameter", "covariate")
)

merged_data <- select(merged_data, gear:parameter_dev, correlation:upper_ci)

plot_multi_filtered <- merged_data

# Plot with reordered levels and filtered data, and adjusted column labels
multi_covariates <- ggplot(
  data = plot_multi_filtered[plot_multi_filtered$covariate != "YEAR", ],
  aes(
    x = covariate_value,
    y = parameter_dev
  )
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(parameter ~ covariate,
    scales = "free",
    labeller = labeller(
      covariate = as_labeller(facet_labels$covariate),
      parameter = as_labeller(facet_labels$parameter)
    )
  ) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "Deviation from global mean") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14)
  ) + # Adjust panel border
  geom_text(
    aes(label = paste(
      round(correlation, 2),
      "(", round(lower_ci, 2), "-",
      round(upper_ci, 2), ")"
    )), # Adjust transparency here
    x = Inf, y = Inf, hjust = 1, vjust = 1, size = 6, color = "darkgray"
  )

print(multi_covariates)


# Do same, but with monofilament data

# Pull in net-level covariates for plotting
covariates <- read.csv("data/plot_covariates.csv")
str(covariates)

covariates$SECCHI <- as.numeric(covariates$SECCHI)

# Center and scale covariates to aid interpreation on plots
covariates$tot_catch_s <- as.numeric(scale(covariates$tot_catch, center = TRUE, scale = TRUE))
covariates$wae_catch_s <- as.numeric(scale(covariates$wae_catch, center = TRUE, scale = TRUE))
covariates$depth_s <- as.numeric(scale(covariates$DEPTH, center = TRUE, scale = TRUE))
covariates$secchi_s <- as.numeric(scale(covariates$SECCHI, center = TRUE, scale = TRUE))

# Model output for bi-normal selectivity function for monofilament data
rep <- readRDS("data/rep_binorm_bysite_mono_all.rds")
obj <- readRDS("data/obj_binorm_bysite_mono_all.rds")
opt <- readRDS("data/opt_binorm_bysite_mono_all.rds")

plot_mono <- covariates[covariates$GEAR == 10, ]

# remove this curve since coordinates unreliable and locations in 2012 where no
# mono nets were fished
plot_mono <- plot_mono[plot_mono$CRNII != 2012005 & !is.na(plot_mono$YEAR), ]

# Parameter deviations for each parameter in selectivity function
plot_mono$k1_dev <- rep$par.random[names(rep$par.random) == "k1_dev"]
plot_mono$k2_dev <- rep$par.random[names(rep$par.random) == "k2_dev"]
plot_mono$k3_dev <- rep$par.random[names(rep$par.random) == "k3_dev"]
plot_mono$k4_dev <- rep$par.random[names(rep$par.random) == "k4_dev"]
plot_mono$c_dev <- rep$par.random[names(rep$par.random) == "c_dev"]

# This is only using standardized covariates
plot_mono_covariates <- plot_mono %>%
  pivot_longer(
    cols = c(tot_catch_s:secchi_s),
    names_to = "covariate",
    values_to = "covariate_value"
  )

plot_mono_pars <- plot_mono %>%
  pivot_longer(
    cols = c(k1_dev, k2_dev, k3_dev, k4_dev, c_dev),
    names_to = "parameter",
    values_to = "parameter_dev"
  )

# Merge the two pivoted data frames
plot_mono_long <- merge(plot_mono_covariates, plot_mono_pars, by = "CRNII")

# Select only columns worth keeping
plot_mono_long <- select(plot_mono_long, CRNII, covariate:covariate_value, parameter:parameter_dev)

plot_mono_long$gear <- "mono"

# Calculate Pearson correlation between standardized covariates
# and parameter deviations.

plot_mono <- plot_mono[complete.cases(plot_mono), ]

# Select columns ending with "_dev" and "_s"
dev_columns <- names(plot_mono)[grep("_dev$", names(plot_mono))]
s_columns <- names(plot_mono)[grep("_s$", names(plot_mono))]

# Initialize an empty data frame to store results
correlation_results <- data.frame()

# Loop through each combination of dev and s columns
for (dev_col in dev_columns) {
  for (s_col in s_columns) {
    # Calculate Pearson correlation
    correlation <- cor(plot_mono[[dev_col]], plot_mono[[s_col]])

    # Calculate Fisher transformation
    fisher_transform <- 0.5 * log((1 + correlation) / (1 - correlation))

    # Compute standard error
    n <- nrow(plot_mono)
    standard_error <- 1 / sqrt(n - 3)

    # Compute z value for 95% confidence interval
    z_value <- qnorm(0.975)

    # Compute margin of error
    margin_of_error <- z_value * standard_error

    # Compute lower and upper bounds of the confidence interval for the Fisher transformation
    lower_bound <- fisher_transform - margin_of_error
    upper_bound <- fisher_transform + margin_of_error

    # Apply inverse Fisher transformation to get the confidence interval for the correlation coefficient
    lower_ci <- (exp(2 * lower_bound) - 1) / (exp(2 * lower_bound) + 1)
    upper_ci <- (exp(2 * upper_bound) - 1) / (exp(2 * upper_bound) + 1)

    # Extract gear and parameter names
    gear <- plot_mono$GEAR[1] # Assuming all rows have the same gear
    parameter <- sub("_dev", "", dev_col)
    covariate <- sub("_s", "", s_col)

    # Create a data frame for the current correlation with confidence intervals
    correlation_df <- data.frame(
      gear = gear,
      parameter = parameter,
      covariate = covariate,
      correlation = correlation,
      lower_ci = lower_ci,
      upper_ci = upper_ci
    )

    # Append to the main results data frame
    correlation_results <- rbind(correlation_results, correlation_df)
  }
}

# Print the results
print(correlation_results)

# For facet labels
facet_labels <- list(
  covariate = c(
    "DEPTH" = "Lake depth (m)",
    "tot_catch" = "Total fishes caught",
    "wae_catch" = "Total walleye caught",
    "SECCHI" = "Secchi depth (m)"
  ),
  parameter = c(
    "k1_dev" = "k1",
    "k2_dev" = "k2",
    "k3_dev" = "k3",
    "k4_dev" = "k4",
    "c_dev" = "c"
  )
)

# If using standardized covariates
facet_labels <- list(
  covariate = c(
    "depth" = "Lake depth (m)",
    "tot_catch" = "Total fishes caught",
    "wae_catch" = "Total walleye caught",
    "secchi" = "Secchi depth (m)"
  ),
  parameter = c(
    "k1_dev" = "k1",
    "k2_dev" = "k2",
    "k3_dev" = "k3",
    "k4_dev" = "k4",
    "c_dev" = "c"
  )
)

# Filter out rows with NA values
plot_mono_filtered <- plot_mono_long[complete.cases(plot_mono_long), ]

# Merge correlation estimates and uncertainty around it to the plotting data

plot_mono_filtered$parameter <- sub("_dev.*", "", plot_mono_filtered$parameter)
plot_mono_filtered$covariate <- sub("_s.*", "", plot_mono_filtered$covariate)

correlation_results$gear <- "mono"

merged_data <- merge(plot_mono_filtered,
  correlation_results,
  by = c("gear", "parameter", "covariate")
)

merged_data <- select(merged_data, gear:parameter_dev, correlation:upper_ci)

plot_mono_filtered <- merged_data

# Plot with reordered levels and filtered data, and adjusted column labels
mono_covariates <- ggplot(
  data = plot_mono_filtered[plot_mono_filtered$covariate != "YEAR", ],
  aes(
    x = covariate_value,
    y = parameter_dev
  )
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(parameter ~ covariate,
    scales = "free",
    labeller = labeller(
      covariate = as_labeller(facet_labels$covariate),
      parameter = as_labeller(facet_labels$parameter)
    )
  ) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "Deviation from global mean") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14)
  ) + # Adjust panel border
  geom_text(
    aes(label = paste(
      round(correlation, 2),
      "(", round(lower_ci, 2), "-",
      round(upper_ci, 2), ")"
    )), # Adjust transparency here
    x = Inf, y = Inf, hjust = 1, vjust = 1, size = 6, color = "darkgray"
  )

print(mono_covariates)
print(multi_covariates)

#### End of script ####
