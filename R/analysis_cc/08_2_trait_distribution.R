library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)
library(vegan)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

# Get species niches
niches <- NULL
for (sp in 1:10) {
  sp_niches <- read_csv(file.path(base_dir, "a2_1", paste0("SpeciesPool", sp, ".csv")),
                        show_col_types = FALSE) %>%
    dplyr::select(-LightResponseA, -LightResponseB, -LightResponseC, -MinWind, -MaxWind, -OptimumWind,
           -DispersalKernelWindEffect) %>%
    mutate(SpeciesPool = sp)

  if (is.null(niches)) {
    niches <- sp_niches
  } else {
    niches <- rbind(niches, sp_niches)
  }
}

# compute sds (NA if column is all NA)
sds <- sapply(niches, sd, na.rm = TRUE)
# columns with zero sd or NA sd
problem_cols <- names(sds)[is.na(sds) | sds == 0]

# PCA
pca_data <- niches %>%
  dplyr::select(-all_of(problem_cols)) %>%
  mutate(
    RangeTemp = MaxTemp - MinTemp,
    RangeHum = MaxHum - MinHum,
    RangeLight = MaxLight - MinLight
  ) %>%
  dplyr::select(-MaxTemp, -MinTemp, -MaxHum, -MinHum, -MaxLight, -MinLight)

# --- Species survival categories

species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

# Get stats
species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario) %>%
  dplyr::select(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year, Position) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = Position
  ) %>%
  mutate(
    diff = `CC` - `No CC`
  )

# Final distribution
species_data_2080_2100 <- species_distr_stats %>%
  filter(Year >= 2080) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    AvgPositionNoCC = mean(`No CC`, na.rm = TRUE),
    AvgPositionCC = mean(`CC`, na.rm = TRUE),
    AvgDiff = mean(diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Shift = ifelse(AvgDiff > 0, "Upward", "Downward")) %>%
  drop_na(.)

# Initial distribution
species_data_1981 <- species_distr_stats %>%
  filter(Year == 1981) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    AvgPositionNoCC = mean(`No CC`, na.rm = TRUE),
    AvgPositionCC = mean(`CC`, na.rm = TRUE),
    AvgDiff = mean(diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Shift = ifelse(AvgDiff > 0, "Upward", "Downward")) %>%
  drop_na(.)

# Join traits with distribution data
species_data_2080_2100 <- left_join(species_data_2080_2100, pca_data, by = c("SpeciesPool", "SpeciesID"))
species_data_2080_2100$community <- "Final"
species_data_1981 <- left_join(species_data_1981, pca_data, by = c("SpeciesPool", "SpeciesID"))
species_data_1981$community <- "Initial"

# Combine both datasets
species_data <- rbind(species_data_1981, species_data_2080_2100)

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/functional_analysis")
dir.create(DirectoryPlots)

# Plot Growth rate vs position of all species alive
pdf(file.path(DirectoryPlots, "growth_rate_vs_position_cc_init_vs_final.pdf"), width = 10, height = 8)
ggplot(species_data, aes(x = GrowthRate, y = AvgPositionCC, color = community)) +
  geom_point() +
  labs(
    x = "Growth Rate",
    y = "Average Position (CC)",
    color = "Community"
  ) +
  theme_minimal()
dev.off()

# - Now only final species

# Filter species that are present in the final community
species_final <- species_data %>%
  filter(community == "Final") %>%
  dplyr::select(SpeciesPool, SpeciesID)

# Keep only species that are in that final community
species_filtered <- species_data %>%
  inner_join(species_final, by = c("SpeciesPool", "SpeciesID"))

# Plot and save as PDF
pdf(file.path(DirectoryPlots, "growth_rate_vs_position_cc_init_vs_final_final_species.pdf"), width = 10, height = 8)
ggplot(species_filtered, aes(x = GrowthRate, y = AvgPositionCC, color = community)) +
  geom_point() +
  labs(
    x = "Growth Rate",
    y = "Average Position (CC)",
    color = "Community"
  ) +
  theme_minimal()
dev.off()

# - Now comparing CC and no CC

pdf(file.path(DirectoryPlots, "growth_rate_vs_position_cc_vs_nocc_final_species.pdf"), width = 10, height = 8)

ggplot(species_filtered, aes(x = GrowthRate, color = community)) +
  geom_point(aes(y = AvgPositionNoCC), alpha = 0.4, size = 2) +  # transparent NoCC
  geom_point(aes(y = AvgPositionCC), alpha = 0.9, size = 2.5) +  # solid CC
  labs(
    x = "Growth Rate",
    y = "Average Position",
    color = "Community"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Growth Rate vs Position (CC and NoCC)")

dev.off()

# Compute correlation between GrowthRate and shift for the final community species
cor_test <- cor.test(species_filtered$GrowthRate[species_filtered$community == "Final"],
                     species_filtered$AvgDiff[species_filtered$community == "Final"],
                     method = "pearson")
cor_test # -> -0.2523081 p-value = 0.0498

# Scatter plot growth rate by initial height for final community species
species_data_2080_2100_init <- species_data_1981 %>%
  dplyr::select(SpeciesPool, SpeciesID, AvgPositionCC, AvgDiff) %>%
  rename(InitialCCPosition = AvgPositionCC, InitialDiff = AvgDiff) %>%
  left_join(species_data_2080_2100, ., by = c("SpeciesPool", "SpeciesID"))

# Scatter plot of growth rate vs shift color coded by initial height for final community species
pdf(file.path(DirectoryPlots, "gr_vs_shift_final_no_cc_pos_color_gr.pdf"), width = 10, height = 8)
ggplot(species_data_2080_2100_init, aes(x = GrowthRate, y = AvgDiff, color = AvgPositionNoCC)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma") +
  labs(
    x = "Growth Rate",
    y = "Shift (CC - Non-CC)",
    color = "Final Position (no CC)"
  ) +
  theme_minimal()
dev.off()

cor.test(species_data_2080_2100_init$InitialCCPosition,
         species_data_2080_2100_init$AvgDiff,method = "pearson") # -> no correlation
cor.test(species_data_2080_2100_init$InitialDiff,
         species_data_2080_2100_init$AvgDiff,method = "pearson") # -> no correlation
cor.test(species_data_2080_2100_init$AvgPositionNoCC,
         species_data_2080_2100_init$AvgDiff,method = "pearson") # -> no correlation

# No please Plot the density of Growth rate for the final community vs initial in species_data
pdf(file.path(DirectoryPlots, "growth_rate_density_cc_init_vs_final.pdf"), width = 10, height = 8)
ggplot(species_data, aes(x = GrowthRate, fill = community)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Growth Rate",
    y = "Density",
    fill = "Community"
  ) +
  theme_minimal()
dev.off()

# - Boxplot of Growth rate by shift category for final community species
pdf(file.path(DirectoryPlots, "growth_rate_boxplot_by_shift_final_species.pdf"), width = 10, height = 8)
ggplot(species_filtered %>%
         filter(community == "Final" & (AvgDiff > 0.5 | AvgDiff < -0.5)),
       aes(x = Shift, y = GrowthRate, fill = Shift)) +
  geom_boxplot() +
  labs(
    x = "Shift Category",
    y = "Growth Rate",
    fill = "Shift"
  ) +
  theme_minimal()
dev.off()