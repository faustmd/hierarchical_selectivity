# Script to produce a two panel figure showing locations of nets in Lake Erie
# and numbers of walleye caught by gear configuration.

library(readxl)

# Load gear comparison data

# Effort and abiotic data
abio <- read_excel("data/gear comparison data 3.0.xlsx",
  sheet = "Abiotic data", na = "."
)

# Remove Partnership (not used), 2010 (nets missing meshes), and
# CRN 2012005 because of suspect coordinates
abio <- abio[abio$GEAR != 15 & abio$YEAR != 2010 & abio$CRNII != 2012005, ]

str(abio)

# Bring in Lake Erie shapefile
library(tidyverse)
library(sf)

lake_erie <- st_read("data/GIS/Lake_Erie_Shoreline.shp")
crop_erie <- st_crop(lake_erie, xmin = -83.5, ymin = 41.339459, xmax = -81, ymax = 42.6)

fig_map_bounds <- st_bbox(c(xmin = -83.8586, ymin = 41.266, xmax = -80.5133, ymax = 42.8241)) %>%
  st_as_sfc() %>% # as a simple feature
  st_set_crs(4326) # define coordinate reference system

# Distill down to only info needed for plotting sites

set_info <- select(abio, CRN:YEAR, GEAR, LAT, LONG, GEAR_FISHED)

set_info %>%
  group_by(CRNII) %>%
  summarize(
    mean_lat = mean(LAT),
    mean_lon = mean(LONG)
  ) %>%
  ungroup() -> set_info2

set_info <- merge(set_info2, set_info, by = c("CRNII"))
set_info <- select(set_info, CRNII:mean_lon, YEAR, GEAR_FISHED)

set_info <- set_info %>%
  st_as_sf(coords = c("mean_lon", "mean_lat"), crs = 4326)

str(set_info)

# Plotting
ggplot(data = crop_erie) +
  geom_sf() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  coord_sf(crs = 4326, xlim = c(-83.5, -81)) + # Adjust xlim as needed
  geom_sf(data = set_info, aes(color = GEAR_FISHED)) +
  scale_x_continuous(breaks = seq(-83.5, -81, by = 1.0)) + # Adjust breaks as needed
  facet_wrap(. ~ YEAR) +
  scale_color_manual(
    name = "Gear fished",
    labels = c(
      "both" = "Both",
      "multi" = "Mulfilament only"
    ),
    values = c("both" = "black", "multi" = "white")
  )
theme(axis.text = element_text(color = "black"))

ggplot(data = crop_erie) +
  geom_sf() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  coord_sf(crs = 4326, xlim = c(-83.5, -81)) + # Adjust xlim as needed
  geom_sf(data = set_info, aes(fill = GEAR_FISHED), shape = 21) + # Modified line
  scale_x_continuous(breaks = seq(-83.5, -81, by = 1)) + # Adjust breaks as needed
  facet_wrap(. ~ YEAR, ncol = 3, nrow = 2) +
  scale_fill_manual(
    name = "Gear fished",
    labels = c(
      "both" = "Both",
      "multi" = "Multifilament only"
    ),
    values = c("both" = "black", "multi" = "white")
  ) +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = "top",
    panel.spacing = unit(2, "lines")
  )


# create figures

bio <- read_excel("data/gear comparison data 3.0.xlsx",
  sheet = "Biological data", na = "."
)

str(bio)

# Remove Partnership (not used), 2010 (nets missing meshes), and
# CRN 2012005 because of suspect coordinates, plus only walleye
wae <- bio[bio$GEAR != 15 &
  bio$YEAR != 2010 &
  bio$CRNII != 2012005 &
  bio$SPECIES == 334, ]

wae <- subset(wae, LENGTH > 0 | !is.na(LENGTH)) # removes some fish missing length data
wae <- subset(wae, MESH < 10) # eliminate nonsense meshes
wae$lbin <- round(wae$LENGTH / 25) * 25 # create 25 mm length bins

# summary stuff for results

catch_by_gear <- wae %>%
  group_by(GEAR) %>%
  summarize(catch = length(GEAR))

catch_by_gear_site <- wae %>%
  group_by(YEAR, CRNII, GEAR) %>%
  summarize(catch = n()) %>%
  ungroup() %>%
  group_by(YEAR, GEAR) %>%
  summarize(mean_catch = mean(catch))

m_len <- wae %>%
  group_by(CRNII, GEAR) %>%
  summarize(mean_length = mean(LENGTH), catch = n()) %>%
  ungroup()

m_len$year <- as.character(m_len$CRNII) %>% substr(1, 4)

ggplot(m_len) +
  geom_bar(aes(x = CRNII, y = catch, fill = as.factor(GEAR)), color = "black", stat = "identity") +
  xlab("Site ID") +
  ylab("Number of walleye caught") +
  scale_fill_manual(
    name = "Gear", labels = c("7" = "Multi", "10" = "Mono"),
    values = c("7" = "gray", "10" = "black")
  ) +
  coord_flip()

p07_bar <- ggplot(m_len[m_len$year == 2007, ]) +
  geom_bar(aes(x = CRNII, y = catch, fill = as.factor(GEAR)), color = "black", stat = "identity") +
  xlab("Site ID") +
  ylab("Number of walleye caught") +
  scale_fill_manual(
    name = "Gear", labels = c("7" = "Multi", "10" = "Mono"),
    values = c("7" = "gray", "10" = "black")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  coord_flip()
p07_bar

p07_map <- ggplot(data = crop_erie) +
  geom_sf() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  coord_sf(crs = 4326) +
  geom_sf(data = set_info[set_info$YEAR == 2007, ]) +
  theme(axis.text = element_text(color = "black"))
p07_map
#### Supplemental figures ####

# Summarize catch at 25-mm length bin for each net
# and gear combination

wae_summ <- wae %>%
  group_by(CRNII, GEAR, lbin) %>%
  summarize(count_by_net = n()) %>%
  ungroup() %>%
  complete(CRNII, GEAR, lbin, fill = list(count_by_net = 0))

wae_summ$year <- as.character(wae_summ$CRNII) %>% substr(1, 4)

ggplot(wae_summ[wae_summ$year == 2007, ]) +
  geom_bar(aes(x = lbin, y = count_by_net, fill = as.factor(GEAR)), color = "black", stat = "identity") +
  facet_wrap(. ~ CRNII) +
  ylab("Number of walleye caught") +
  xlab("Length bin (25 mm)") +
  theme_bw() +
  scale_fill_manual(
    name = "Gear", labels = c("7" = "Multifilament", "10" = "Monofilament"),
    values = c("7" = "gray", "10" = "black")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top",
    legend.title = element_text(color = "black", size = 18),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(
      t = 0.5,
      r = 0.5,
      b = 0.5,
      l = 0.5,
      unit = "in"
    )
  ) # Adjust panel border)

ggplot(wae_summ[wae_summ$year == 2008, ]) +
  geom_bar(aes(x = lbin, y = count_by_net, fill = as.factor(GEAR)), color = "black", stat = "identity") +
  facet_wrap(. ~ CRNII) +
  ylab("Number of walleye caught") +
  xlab("Length bin (25 mm)") +
  theme_bw() +
  scale_fill_manual(
    name = "Gear", labels = c("7" = "Multifilament", "10" = "Monofilament"),
    values = c("7" = "gray", "10" = "black")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top",
    legend.title = element_text(color = "black", size = 18),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(
      t = 0.5,
      r = 0.5,
      b = 0.5,
      l = 0.5,
      unit = "in"
    )
  ) # Adjust panel border)

ggplot(wae_summ[wae_summ$year == 2011, ]) +
  geom_bar(aes(x = lbin, y = count_by_net, fill = as.factor(GEAR)), color = "black", stat = "identity") +
  facet_wrap(. ~ CRNII) +
  ylab("Number of walleye caught") +
  xlab("Length bin (25 mm)") +
  theme_bw() +
  scale_fill_manual(
    name = "Gear", labels = c("7" = "Multifilament", "10" = "Monofilament"),
    values = c("7" = "gray", "10" = "black")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top",
    legend.title = element_text(color = "black", size = 18),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(
      t = 0.5,
      r = 0.5,
      b = 0.5,
      l = 0.5,
      unit = "in"
    )
  ) # Adjust panel border)

ggplot(wae_summ[wae_summ$year == 2012, ]) +
  geom_bar(aes(x = lbin, y = count_by_net, fill = as.factor(GEAR)), color = "black", stat = "identity") +
  facet_wrap(. ~ CRNII) +
  ylab("Number of walleye caught") +
  xlab("Length bin (25 mm)") +
  theme_bw() +
  scale_fill_manual(
    name = "Gear", labels = c("7" = "Multifilament", "10" = "Monofilament"),
    values = c("7" = "gray", "10" = "black")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top",
    legend.title = element_text(color = "black", size = 18),
    panel.spacing = unit(1.5, "lines"),
    plot.margin = margin(
      t = 0.5,
      r = 0.5,
      b = 0.5,
      l = 0.5,
      unit = "in"
    )
  ) # Adjust panel border)

ggplot(wae_summ[wae_summ$year == 2013, ]) +
  geom_bar(aes(x = lbin, y = count_by_net, fill = as.factor(GEAR)), color = "black", stat = "identity") +
  facet_wrap(. ~ CRNII) +
  ylab("Number of walleye caught") +
  xlab("Length bin (25 mm)") +
  theme_bw() +
  scale_fill_manual(
    name = "Gear", labels = c("7" = "Multifilament", "10" = "Monofilament"),
    values = c("7" = "gray", "10" = "black")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", size = 1), # Adjust strip background
    panel.border = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 14),
    legend.text = element_text(color = "black", size = 18),
    legend.position = "top",
    legend.title = element_text(color = "black", size = 18),
    panel.spacing = unit(1.5, "lines"),
    plot.margin = margin(
      t = 0.5,
      r = 0.5,
      b = 0.5,
      l = 0.5,
      unit = "in"
    )
  ) # Adjust panel border)

#### End of script ####
