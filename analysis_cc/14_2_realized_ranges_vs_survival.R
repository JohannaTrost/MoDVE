library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)
library(vegan)
library(tibble)
library(randomForest)
library(patchwork)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

# --- Get species stats
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

# --- Survival

species_survival_cat <- species_distr_stats %>%
  filter(Year >= 2080) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    has_NoCC = any(!is.na(`No CC`)),
    has_CC   = any(!is.na(CC)),
    .groups = "drop"
  ) %>%
  mutate(
    survival = case_when(
      has_CC & !has_NoCC ~ "Survived CC",
      !has_CC & has_NoCC ~ "Died with CC",
      has_CC & has_NoCC  ~ "Survived CC",
      TRUE ~ "none"  # optional: if both are NA
    )
  ) %>%
  dplyr::select(SpeciesPool, SpeciesID, survival)

# --- Extract env. conditions

# - 1. Get time and location of individuals od interest
filtered_species_distr <- species_distr_stats %>%
  group_by(SpeciesPool, SpeciesID) %>%
  arrange(Year, .by_group = TRUE) %>%  # ensure data is sorted by Year
  reframe(
    diff = mean(tail(diff[!is.na(diff)], 20), na.rm = TRUE),
    CC = mean(tail(CC[!is.na(CC)], 20), na.rm = TRUE),
    last_year_alive = max(Year[!is.na(diff)]),
  ) %>%
  #filter(last_year_alive >= 2050) %>%
  select(SpeciesPool, SpeciesID, last_year_alive, diff) %>%
  right_join(species_distr, ., by = c("SpeciesPool", "SpeciesID")) %>%
  filter(Year <= last_year_alive & Year >= last_year_alive - 20) %>%
  drop_na(.) # Drop all species that only survived in the sim. with community

# - 2. Loop over years and extract MC matrices and env. conditions
scenario_dir <- list("CC" = "climdata_era5_cmip6_1981-2100_ssp245",
                     "No CC" = "climdata_era5_cmip6_1981-2100_ssp245_no_cc")
species_env <- NA
for (scenario in c("CC", "No CC")) {
  for (forest in 0:2) {
    species_forest <- filtered_species_distr %>% filter(ForestID == forest, Scenario == scenario)

    for (ts in sort(unique(species_forest$TimeStep))) {
      species_locations <- species_forest %>% filter(TimeStep == ts)

      # Load mC matrix
      mh <- readRDS(file.path(base_dir, scenario_dir[[scenario]], "a1_2",
                              paste0("forest", forest), paste0("microhabitatMatrix", ts, ".rds")))
      mc <- mh[,,, c(3,4,5)] # Get microclimate
      substrate_area <- mh[,,, 1] # Get substrate area

      # For each individual, get the max non-zero value in substrate_area[X, Y, ]
      species_locations$VegetationHeight <- mapply(function(x, y) {
          layer_values <- substrate_area[x, y, ]
          max_non_zero <- max(which(layer_values != 0))
        }, species_locations$X, species_locations$Y)

      # Get the env conditions at the species position
      species_locations[, c("Light", "Hum", "Temp")] <- t(mapply(function(x, y, z) {
        mc[x, y, z, ]
      }, species_locations$X, species_locations$Y, species_locations$Z))

      # Get the minimum potential env. conditions on the vertical axes
      species_locations[, c("MinPotLight", "MinPotHum", "MinPotTemp")] <- t(
        mapply(function(x, y, hgt) {
          if (hgt <= 1) {
            mc[x, y, 1, ]  # Return the values at height 1 if hgt <= 1
          } else {
            apply(mc[x, y, 1:hgt, ], 2, min)  # Otherwise, return the min values across 1:hgt
          }
        }, species_locations$X, species_locations$Y, species_locations$VegetationHeight)
      )

      # Get the max potential env. conditions on the vertical axes
      species_locations[, c("MaxPotLight", "MaxPotHum", "MaxPotTemp")] <- t(
        mapply(function(x, y, hgt) {
          if (hgt <= 1) {
            mc[x, y, 1, ]  # Return the values at height 1 if hgt <= 1
          } else {
            apply(mc[x, y, 1:hgt, ], 2, max)  # Otherwise, return the min values across 1:hgt
          }
        }, species_locations$X, species_locations$Y, species_locations$VegetationHeight)
      )

      # Specifiy scenario
      species_locations$Scenario <- scenario

      if (all(is.na(species_env))) {
        species_env <- species_locations
      } else {
        species_env <- rbind(species_locations, species_env)
      }
    }
  }
}

# Now get the theoretical niche specifications
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

# Species env. conditions
species_env_agg <- species_env %>%
  group_by(Scenario, SpeciesPool, SpeciesID) %>%
  summarize(AvgTempCond = mean(Temp, na.rm = TRUE),
            MaxTempCond = max(Temp, na.rm = TRUE),
            MinTempCond = min(Temp, na.rm = TRUE),
            AvgHumCond = mean(Hum, na.rm = TRUE),
            MaxHumCond = max(Hum, na.rm = TRUE),
            MinHumCond = min(Hum, na.rm = TRUE),
            AvgLightCond = mean(Light, na.rm = TRUE) * 900,
            MaxLightCond = max(Light, na.rm = TRUE) * 900,
            MinLightCond = min(Light, na.rm = TRUE) * 900,
            # Potential env. conditions per species
            MaxPotLight = max(MaxPotLight, na.rm = TRUE) * 900,
            MaxPotHum = max(MaxPotHum, na.rm = TRUE),
            MaxPotTemp = max(MaxPotTemp, na.rm = TRUE),
            MinPotLight = min(MinPotLight, na.rm = TRUE) * 900,
            MinPotHum = min(MinPotHum, na.rm = TRUE),
            MinPotTemp = min(MinPotTemp, na.rm = TRUE),
            AvgDiff = mean(diff, na.rm = TRUE),
            AvgHeight = mean(VegetationHeight, na.rm = TRUE),
            .groups = "drop")

# Merge realized and fundamental niches
species_env_niches <- left_join(species_env_agg, niches, by = c("SpeciesID", "SpeciesPool"))

# --- Constrain realizable values by potential niche

species_env_niches_c <- species_env_niches %>% # Constrain potential limits with fundamental niche
  mutate(
    MaxSuitTemp = pmin(MaxPotTemp, MaxTemp),
    MinSuitTemp = pmax(MinPotTemp, MinTemp),
    MaxSuitHum = pmin(MaxPotHum, MaxHum),
    MinSuitHum = pmax(MinPotHum, MinHum),
    MaxSuitLight = pmin(MaxPotLight, MaxLight),
    MinSuitLight = pmax(MinPotLight, MinLight),
  )

# --- Compute range sizes

range_sizes <- species_env_niches_c %>%
  mutate(TempRealizedSize = MaxTempCond - MinTempCond,
         HumRealizedSize = MaxHumCond - MinHumCond,
         LightRealizedSize = MaxLightCond - MinLightCond,
         TempAvailRangeSize = MaxPotTemp - MinPotTemp,
         HumAvailRangeSize = MaxPotHum - MinPotHum,
         LightAvailRangeSize = MaxPotLight - MinPotLight,
         TempRangeSize = MaxSuitTemp - MinSuitTemp,
         HumRangeSize = MaxSuitHum - MinSuitHum,
         LightRangeSize = MaxSuitLight - MinSuitLight,
         TempNicheSize = MaxTemp - MinTemp,
         HumNicheSize = MaxHum - MinHum,
         LightNicheSize = MaxLight - MinLight) %>%
  select(Scenario, SpeciesPool, SpeciesID, TempRealizedSize, HumRealizedSize, LightRealizedSize, TempRangeSize,
         HumRangeSize, LightRangeSize, TempNicheSize, HumNicheSize, LightNicheSize, TempAvailRangeSize,
         HumAvailRangeSize, LightAvailRangeSize
  )

# --- Merge survival with ranges/niches data

range_sizes_survival <- left_join(species_survival_cat, range_sizes, by = c("SpeciesID", "SpeciesPool"))

# Scatterplot
DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc")
pdf(file.path(DirectoryPlots, "fundamental_vs_realized_by_survival.pdf"))
# Calculate axis limits (same for x and y)
axis_limits <- c(0, max(range_sizes_survival$HumNicheSize, range_sizes_survival$HumRealizedSize, na.rm = TRUE))

ggplot(data = range_sizes_survival %>% filter(Scenario == "CC"),
       aes(x = HumNicheSize, y = HumRealizedSize, color = survival)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # 1:1 line
  coord_fixed() +  # Forces equal scaling for x and y
  xlim(axis_limits) +  # Set equal limits
  ylim(axis_limits) +
  labs(x = "Humidity fundamental niche size", y = "Humidity realized size",
       color = "Survival Status") +
  theme_minimal() +
  scale_color_manual(values = c("red", "green"))
dev.off()

# - Show PDFs

# Pivot the data to long format for faceting
plot_data <- range_sizes_survival %>%
  filter(Scenario == "No CC") %>%
  pivot_longer(
    cols = c(HumRealizedSize, HumRangeSize, HumNicheSize),
    names_to = "Metric",
    values_to = "Value"
  )

# Calculate global axis limits (same for all subplots)
global_limits <- c(
  0,
  max(plot_data$Value, na.rm = TRUE) * 1.05  # Add 5% padding
)

pdf(file.path(DirectoryPlots, "ranges_niches_distr_by_survival_baseline.pdf"))
# Create the faceted density plot
ggplot(plot_data, aes(x = Value, fill = survival)) +
  geom_density(alpha = 0.5, adjust = 1.5) +  # Adjust bandwidth as needed
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +  # One subplot per metric
  coord_cartesian(xlim = global_limits) +  # Equal x-axis limits
  labs(x = "Value", y = "Density", fill = "Survival Status") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green")) +
  theme(legend.position = "top")  # Move legend to top
dev.off()

# --- Compare available ranges -> for this uncomment the filter for all > 2050

# Pivot the data to long format for faceting
plot_data <- range_sizes %>%
  pivot_longer(
    cols = c(TempRealizedSize, TempRangeSize, TempNicheSize, TempAvailRangeSize),
    names_to = "Metric",
    values_to = "Value"
  )

# Calculate global axis limits (same for all subplots)
global_limits <- c(
  0,
  max(plot_data$Value, na.rm = TRUE) * 1.05  # Add 5% padding
)

pdf(file.path(DirectoryPlots, "ranges_niches_distr_by_scenario.pdf"))
# Create the faceted density plot
ggplot(plot_data, aes(x = Value, fill = Scenario)) +
  geom_density(alpha = 0.5, adjust = 1.5) +  # Adjust bandwidth as needed
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +  # One subplot per metric
  coord_cartesian(xlim = global_limits) +  # Equal x-axis limits
  labs(x = "Value", y = "Density", fill = "Scenario") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue")) +
  theme(legend.position = "top")  # Move legend to top
dev.off()