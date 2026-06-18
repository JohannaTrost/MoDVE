library(readr)
library(ggplot2)
library(dplyr)
library(patchwork)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")
DirectoryPlots <- "../../figs/a5_plots_test/cc_vs_no_cc/position_shift"
indiv_env_niches <- read_csv(file.path(DirectoryPlots, "a5_individuals_cc_vs_no_cc_niches_envCond_shift_v2.csv"))

# - Compute distance to tolerances
distances <- indiv_env_niches %>%
  mutate(MaxTempDist = abs(MaxTemp - Temp),
         MaxHumDist = abs(MaxHum - Hum),
         MinLightDist = abs(MinLight - Light),
         OptTempDist = abs(OptimumTemp - Temp),
         OptHumDist = abs(OptimumHum - Hum),
         OptLightDist = abs(OptimumLight - Light)
  ) %>%
  group_by(Scenario, SpeciesPool, SpeciesID) %>%
  summarise(MaxTempDist = mean(MaxTempDist, na.rm = TRUE),
            MaxHumDist = mean(MaxHumDist, na.rm = TRUE),
            MinLightDist = mean(MinLightDist, na.rm = TRUE),
            OptTempDist = mean(OptTempDist, na.rm = TRUE),
            OptHumDist = mean(OptHumDist, na.rm = TRUE),
            OptLightDist = mean(OptLightDist, na.rm = TRUE),
  )

# - Distance to max. limit

# Reshape data to long format for ggplot
max_dists_long <- distances %>%
  select(Scenario, SpeciesPool, SpeciesID, MaxTempDist, MaxHumDist, MinLightDist) %>%
  pivot_longer(
    cols = c(MaxTempDist, MaxHumDist, MinLightDist),
    names_to = "Variable",
    values_to = "DistanceToMax"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      MaxTempDist = "Temperature",
                      MaxHumDist = "Relative humidity",
                      MinLightDist = "Light")
  )

# - Distance to optimum max. limit

# Reshape data to long format for ggplot
opt_dists_long <- distances %>%
  select(Scenario, SpeciesPool, SpeciesID, OptTempDist, OptHumDist, OptLightDist) %>%
  pivot_longer(
    cols = c(OptTempDist, OptHumDist, OptLightDist),
    names_to = "Variable",
    values_to = "DistanceToOpt"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      OptTempDist = "Temperature",
                      OptHumDist = "Relative humidity",
                      OptLightDist = "Light")
  )

# Merge and save distance df
dist_df <- inner_join(opt_dists_long, max_dists_long, by = c("Scenario", "SpeciesPool", "SpeciesID", "Variable"))
write_csv(dist_df, file.path(base_dir, "a5_opt_max_distances_cc_vs_no_cc_long_v1.csv"))

# Define colors for scenarios
colors <- c('CC' = '#f7766e', 'No CC' = '#004aad')

niche_filling <- read_csv(file.path(base_dir, "a5_niche_range_filling_cc_vs_no_cc_v3.csv"))

# Create the boxplot
output_file2 <- file.path(DirectoryPlots, "boxplot_niche_filling.pdf")
#pdf(output_file2, width = 7, height = 6)
ggplot(niche_filling, aes(x = Variable, y = Value, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  labs(
    x = "",
    y = "Niche filling (%)",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
#dev.off()

# --- Combined plots

colors <- c(
  "Baseline" = "#004aad",
  "Climate change" = "#f7766e"
)

make_boxplot <- function(df, yvar, ylab, legend = "none", ylim = NULL) {

  df$Scenario <- factor(
    df$Scenario,
    levels = c("No CC", "CC"),
    labels = c("Baseline", "Climate change")
  )

  p <- ggplot(df, aes(x = Scenario, y = .data[[yvar]], fill = Scenario)) +
    geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
    scale_fill_manual(values = colors, name = NULL) +
    labs(x = "", y = ylab) +
    theme_minimal() +
    theme(
      text = element_text(size = 26),
      axis.text.y = element_text(size = 23),
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = legend,
      plot.title = element_blank()  # Ensures no plot title
    )

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }

  return(p)
}

# ---- 1. Compute global ranges for shared scales ----

range_T <- range(
  c(max_dists_long$DistanceToMax[max_dists_long$Variable == "Temperature"],
    opt_dists_long$DistanceToOpt[opt_dists_long$Variable == "Temperature"]),
  na.rm = TRUE
)
range_RH <- range(
  c(max_dists_long$DistanceToMax[max_dists_long$Variable == "Relative humidity"],
    opt_dists_long$DistanceToOpt[opt_dists_long$Variable == "Relative humidity"]),
  na.rm = TRUE
)
range_L <- range(
  c(max_dists_long$DistanceToMax[max_dists_long$Variable == "Light"],
    opt_dists_long$DistanceToOpt[opt_dists_long$Variable == "Light"]),
  na.rm = TRUE
)

# ---- 2. Build plots (add shared scales where needed) ----

# Upper limit
p_max_T <- make_boxplot(
  filter(max_dists_long, Variable == "Temperature"),
  "DistanceToMax", "Distance to\nmax. temperature (°C)"
) + coord_cartesian(ylim = range_T)

p_max_RH <- make_boxplot(
  filter(max_dists_long, Variable == "Relative humidity"),
  "DistanceToMax", "Distance to\nmax. humidity (%)"
) + coord_cartesian(ylim = range_RH)

p_max_Light <- make_boxplot(
  filter(max_dists_long, Variable == "Light"),
  "DistanceToMax",
  # multiline label: top line "Distance to", second line with unit including μmol m^-2 s^-1
  expression(atop("Distance to min.",
                  paste("light (", mu, "mol ", m^{-2}, " ", s^{-1}, ")"))),
  ylim = range_L
)
# Optimum
p_opt_T <- make_boxplot(
  filter(opt_dists_long, Variable == "Temperature"),
  "DistanceToOpt", "Distance to\ntemperature opt. (°C)"
) + coord_cartesian(ylim = range_T)

p_opt_RH <- make_boxplot(
  filter(opt_dists_long, Variable == "Relative humidity"),
  "DistanceToOpt", "Distance to\nhumidity opt. (%)"
) + coord_cartesian(ylim = range_RH)

p_opt_Light <- make_boxplot(
  filter(opt_dists_long, Variable == "Light"), "DistanceToOpt",
  expression(atop("Distance to light",
                  paste("opt. (", mu, "mol ", m^{-2}, " ", s^{-1}, ")")))
) + coord_cartesian(ylim = range_L)

# FN filling (percentage units differ → separate scales)
p_nf_T <- make_boxplot(
  filter(niche_filling, Variable == "Temperature FN filling"),
  "Value", "Temperature\nFN filling (%)", ylim = c(0, 100)
)
p_nf_RH <- make_boxplot(
  filter(niche_filling, Variable == "Relative humidity FN filling"),
  "Value", "Humidity\nFN filling (%)", ylim = c(0, 100)
)
p_nf_Light <- make_boxplot(
  filter(niche_filling, Variable == "Light FN filling"),
  "Value", "Light\nFN filling (%)", ylim = c(0, 100),
  legend = TRUE
)

# PN filling
p_rf_T <- make_boxplot(
  filter(niche_filling, Variable == "Temperature PN filling"),
  "Value", "Temperature\nPN filling (%)", ylim = c(0, 100)
)
p_rf_RH <- make_boxplot(
  filter(niche_filling, Variable == "Relative humidity PN filling"),
  "Value", "Humidity\nPN filling (%)", ylim = c(0, 100)
)
p_rf_Light <- make_boxplot(
  filter(niche_filling, Variable == "Light PN filling"),
  "Value", "Light\nPN filling (%)", ylim = c(0, 100),
  legend = TRUE
)

# PN size
p_pns_T <- make_boxplot(
  filter(niche_filling, Variable == "Temperature PN size"),
  "Value", "Temperature\nPN size (°C)"
)
p_pns_RH <- make_boxplot(
  filter(niche_filling, Variable == "Relative humidity PN size"),
  "Value", "Humidity\nPN size (%)"
)
p_pns_Light <- make_boxplot(
  filter(niche_filling, Variable == "Light PN size"),
  "Value", "Light\nPN size (lux)",
  legend = TRUE
)

# Geographic range filling
p_vrf <- make_boxplot(
  filter(niche_filling, Variable == "Vertical range filling"),
  "Value", "Vertical range\nfilling (%)"
)
p_hrf <- make_boxplot(
  filter(niche_filling, Variable == "Horizontal range filling"),
  "Value", "Horizontal\nrange filling (%)"
)
p_3drf <- make_boxplot(
  filter(niche_filling, Variable == "3D range filling"),
  "Value", "3D range\nfilling (%)",
  legend = TRUE
)

combined_plot <-
  ( (p_vrf | p_hrf | p_3drf) /
    (p_nf_T | p_rf_T | p_pns_T) /
    (p_nf_RH | p_rf_RH | p_pns_RH) /
    (p_nf_Light | p_rf_Light | p_pns_Light) /
    (p_max_T | p_max_RH | p_opt_T)
  ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.key.width = unit(1.5, "cm")
  )


# ---- 4. Save ----
output_file <- file.path(DirectoryPlots, "niche_range_filling_v4.pdf")
ggsave(output_file, combined_plot, width = 14, height = 16)

## --- try sth

niche_filling$ScenarioNames <- factor(
  niche_filling$Scenario,
  levels = c("No CC", "CC"),
  labels = c("Baseline", "Climate change")
)

# save as csv
write_csv(niche_filling, file.path(base_dir, "a5_niche_range_filling_cc_vs_no_cc_long_v3.csv"))

# Prepare FN/PN filling data
fn_pn_data <- niche_filling %>%
  filter(grepl("FN filling|PN filling", Variable)) %>%
  mutate(
    FillingType = ifelse(grepl("FN", Variable), "FN", "PN"),
    Variable = gsub(" FN filling| PN filling", "", Variable),
    Group = paste0(Variable, " ", FillingType, " filling")
  )

# Factor the Group column to ensure correct ordering
fn_pn_data$Group <- factor(
  fn_pn_data$Group,
  levels = c(
    "Temperature FN filling", "Temperature PN filling",
    "Relative humidity FN filling", "Relative humidity PN filling",
    "Light FN filling", "Light PN filling"
  )
)

# Plot geo-range filling
p_geo_rf <- niche_filling %>%
  filter(Variable == "Vertical range filling" | Variable == "Horizontal range filling" |
         Variable == "3D range filling") %>%
  ggplot(., aes(x = Variable, y = Value, fill = ScenarioNames)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.alpha = 0.5,
    notch = TRUE,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(
    values = colors,
    name = NULL
  ) +
  labs(
    x = "",
    y = "Range filling (%)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 26),
    axis.text.y = element_text(size = 23),
    axis.text.x = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 20)
  )

# Plot PN and FN filling for each env. axis
p_fn_pn <- ggplot(fn_pn_data, aes(x = Group, y = Value, fill = ScenarioNames)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.alpha = 0.5,
    notch = TRUE,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(
    values = colors,
    name = NULL
  ) +
  labs(
    x = "",
    y = "Niche filling (%)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 26),
    axis.text.y = element_text(size = 23),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 20)
  )

final_plot <- (
  # Row 1: Geographical Range Filling
  (p_geo_rf + p_max_T) +
  # Row 2: FN/PN Filling
  (p_fn_pn + p_max_RH) +
  # Row 3: PN Size + Distance to Optimum
  (p_pns_T + p_pns_RH + p_pns_Light + p_opt_T)
) +
  plot_layout(
    guides = "collect",
    ncol = 4,
    widths = c(2, 2, 1, 1),  # p_geo_rf and p_fn_pn span 2 columns
    heights = c(1, 1, 1)
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.key.width = unit(1.5, "cm")
  )

output_file2 <- file.path(DirectoryPlots, "niche_range_filling_v5.pdf")
ggsave(output_file2, final_plot, width = 12, height = 20)
