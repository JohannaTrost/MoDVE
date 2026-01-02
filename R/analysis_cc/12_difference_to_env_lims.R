library(readr)

DirectoryPlots <- "../../figs/a5_plots_test/cc_vs_no_cc/position_shift"
indiv_env_niches <- read_csv(file.path(DirectoryPlots, "a5_individuals_cc_vs_no_cc_niches_envCond_shift.csv"))

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

# - Plot distance to max max. limit

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

# Define colors for scenarios
colors <- c('CC' = '#f7766e', 'No CC' = '#004aad')

# Create the boxplot
output_file2 <- file.path(DirectoryPlots, "density_absdistance_to_max_tol.pdf")
pdf(output_file2, width = 7, height = 9)  # taller height for 3 rows

ggplot(max_dists_long, aes(x = DistanceToMax, fill = Scenario)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ Variable, ncol = 1, scales = "free") +  # 3 rows, one per variable
  labs(
    x = "Distance to max. limit",
    y = "Density",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 14)
  )

dev.off()

# - Plot distance to optimum max. limit

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

# Create the boxplot
output_file2 <- file.path(DirectoryPlots, "density_absdistance_to_opt_tol.pdf")
pdf(output_file2, width = 7, height = 9)  # taller height for 3 rows

ggplot(opt_dists_long, aes(x = DistanceToOpt, fill = Scenario)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ Variable, ncol = 1, scales = "free") +  # 3 rows, one per variable
  labs(
    x = "Distance to optimum value",
    y = "Density",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 14)
  )

dev.off()

# ---- boxplots

output_file2 <- file.path(DirectoryPlots, "boxplot_distance_to_max_tol_3col.pdf")
pdf(output_file2, width = 12, height = 5)  # wider for 3 columns

ggplot(max_dists_long, aes(x = Scenario, y = DistanceToMax, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ Variable, ncol = 3, scales = "free_y") +  # 3 columns, one per variable
  labs(
    x = "",
    y = "Distance to max. limit",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14)
  )

dev.off()

# - Opt
output_file2 <- file.path(DirectoryPlots, "boxplot_distance_to_opt_tol_3col.pdf")
pdf(output_file2, width = 12, height = 5)  # wider for 3 columns

ggplot(opt_dists_long, aes(x = Scenario, y = DistanceToOpt, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ Variable, ncol = 3, scales = "free_y") +  # 3 columns, one per variable
  labs(
    x = "",
    y = "Distance to the species optimum",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14)
  )

dev.off()


# Define colors for scenarios
colors <- c('CC' = '#f7766e', 'No CC' = '#004aad')

niche_filling <- read_csv(file.path(base_dir, "a5_niche_filling_cc_vs_no_cc.csv"))

# Create the boxplot
output_file2 <- file.path(DirectoryPlots, "boxplot_niche_filling.pdf")
pdf(output_file2, width = 7, height = 6)
ggplot(niche_filling, aes(x = Variable, y = NicheFill, fill = Scenario)) +
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
dev.off()

# --- Combined plots

colors <- c(
  "Baseline" = "#004aad",
  "Climate change" = "#f7766e"
)


make_boxplot <- function(df, yvar, ylab, legend = "none") {

  df$Scenario <- factor(
    df$Scenario,
    levels = c("No CC", "CC"),
    labels = c("Baseline", "Climate change")
  )
  ggplot(df, aes(x = Scenario, y = .data[[yvar]], fill = Scenario)) +
    geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
    scale_fill_manual(values = colors, name = NULL) +
    facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
    labs(x = "", y = ylab) +
    theme_minimal() +
    theme(
      text = element_text(size = 26),
      axis.text.y = element_text(size = 23),
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = legend
    )
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
                  paste("light (", mu, "mol ", m^{-2}, " ", s^{-1}, ")")))
) + coord_cartesian(ylim = range_L)
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
  filter(niche_filling, Variable == "Temperature"),
  "NicheFill", "Temperature\nniche filling (%)"
)
p_nf_RH <- make_boxplot(
  filter(niche_filling, Variable == "Relative humidity"),
  "NicheFill", "Humidity\nniche filling (%)"
)
p_nf_Light <- make_boxplot(
  filter(niche_filling, Variable == "Light"),
  "NicheFill", "Light\nniche filling (%)",
  legend = TRUE
)

combined_plot <-
  (
    (p_max_T | p_opt_T | p_nf_T) /
    (p_max_RH | p_opt_RH | p_nf_RH) /
    (p_max_Light | p_opt_Light | p_nf_Light)
  ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.key.width = unit(1.5, "cm")
  )


# ---- 4. Save ----
output_file <- file.path(DirectoryPlots, "niche_filling_v2.pdf")
ggsave(output_file, combined_plot, width = 14, height = 14)

