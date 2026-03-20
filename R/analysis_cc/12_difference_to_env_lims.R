library(readr)

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

# Define colors for scenarios
colors <- c('CC' = '#f7766e', 'No CC' = '#004aad')

niche_filling <- read_csv(file.path(base_dir, "a5_niche_range_filling_cc_vs_no_cc_v2.csv"))

# Create the boxplot
output_file2 <- file.path(DirectoryPlots, "boxplot_niche_filling.pdf")
#pdf(output_file2, width = 7, height = 6)
ggplot(niche_filling, aes(x = Variable, y = Fill, fill = Scenario)) +
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

# Niche filling (percentage units differ → separate scales)
p_nf_T <- make_boxplot(
  filter(niche_filling, Variable == "Temperature niche fill"),
  "Fill", "Temperature\nniche filling (%)", ylim = c(0, 100)
)
p_nf_RH <- make_boxplot(
  filter(niche_filling, Variable == "Relative humidity niche fill"),
  "Fill", "Humidity\nniche filling (%)", ylim = c(0, 100)
)
p_nf_Light <- make_boxplot(
  filter(niche_filling, Variable == "Light niche fill"),
  "Fill", "Light\nniche filling (%)", ylim = c(0, 100),
  legend = TRUE
)

# Range filling
p_rf_T <- make_boxplot(
  filter(niche_filling, Variable == "Temperature range fill"),
  "Fill", "Temperature\nrange filling (%)", ylim = c(0, 100)
)
p_rf_RH <- make_boxplot(
  filter(niche_filling, Variable == "Relative humidity range fill"),
  "Fill", "Humidity\nrange filling (%)", ylim = c(0, 100)
)
p_rf_Light <- make_boxplot(
  filter(niche_filling, Variable == "Light range fill"),
  "Fill", "Light\nrange filling (%)", ylim = c(0, 100),
  legend = TRUE
)

combined_plot <-
  (
    (p_max_T | p_nf_T | p_rf_T) /
    (p_max_RH | p_nf_RH | p_rf_RH) /
    (p_opt_T | p_nf_Light | p_rf_Light)
  ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.key.width = unit(1.5, "cm")
  )


# ---- 4. Save ----
output_file <- file.path(DirectoryPlots, "niche_range_filling_v3.pdf")
ggsave(output_file, combined_plot, width = 14, height = 14)

