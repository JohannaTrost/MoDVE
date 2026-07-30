# -----
# 1. Compute fundamental, potential, realized and geographical niches
# 2. Relate niche shifts to position shifts
# 3. Visual comparison of range sizes

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(patchwork) # for combining plots
library(purrr)

# Define a function to compute results for a single species
compute_suitable_spaces <- function(mc, substrate_area, currSpeciesTraits, i) {
  minL <- currSpeciesTraits[i,]$MinLightScaled
  maxL <- currSpeciesTraits[i,]$MaxLightScaled
  minH <- currSpeciesTraits[i,]$MinHum
  maxH <- currSpeciesTraits[i,]$MaxHum
  minT <- currSpeciesTraits[i,]$MinTemp
  maxT <- currSpeciesTraits[i,]$MaxTemp
  MeanAreaOccupied <- currSpeciesTraits[i,]$MeanAreaOccupied

  lightMask <- mc[,,, 1] >= minL & mc[,,, 1] <= maxL
  humMask <- mc[,,, 2] >= minH & mc[,,, 2] <= maxH
  tempMask <- mc[,,, 3] >= minT & mc[,,, 3] <= maxT
  areaMask <- substrate_area >= MeanAreaOccupied

  suitableMask <- (lightMask & humMask & tempMask & areaMask)

  suitableHeights <- apply(suitableMask, 3, sum)
  suitableHeights[is.na(suitableHeights)] <- 0
  suitableColumns <- apply(suitableMask, c(1,2), any)  # TRUE if ANY height suitable at (X,Y)

  list(
    nSuitableHeights = sum(suitableHeights != 0),
    nSuitableVoxels = sum(suitableMask, na.rm = TRUE),
    nSuitableLight = sum(lightMask, na.rm = TRUE),
    nSuitableHum = sum(humMask, na.rm = TRUE),
    nSuitableTemp = sum(tempMask, na.rm = TRUE),
    nSuitableArea = sum(areaMask, na.rm = TRUE),
    nSuitableColumns = sum(suitableColumns, na.rm = TRUE)
  )
}

# Directories
base_dir <- file.path("../modve_data/modve_output/regua")
DirectoryPlots <- file.path("../modve_figs/climate_change/niches")

dir.create(DirectoryPlots, recursive = TRUE)

# Load species distribution
species_distr <- read_csv(file.path(base_dir, "species_distribution_cc_vs_no_cc.csv"))

# Get position, range and species shift (and widen)
species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
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
  ) %>%
  drop_na(.)

# Summarize shift and position under CC across forests and years
summary_shift <- species_distr_stats %>%
  group_by(SpeciesPool, SpeciesID) %>%
  arrange(Year, .by_group = TRUE) %>%  # ensure data is sorted by Year
  reframe(
    diff = mean(tail(diff[!is.na(diff)], 20), na.rm = TRUE),
    CC = mean(tail(CC[!is.na(CC)], 20), na.rm = TRUE),
    last_year_alive = max(Year[!is.na(diff)]),
  ) %>%
  filter(last_year_alive >= 2050) %>%
  mutate(Shift = ifelse(diff > 0, "Upward", "Downward"))

# ---- Get the theoretical niche specifications

niches <- NULL
for (sp in 1:10) {
  sp_niches <- read_csv(file.path(base_dir, "species_pools", paste0("SpeciesPool", sp, ".csv")),
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

# ---- Extract env. conditions

SurfaceBiomassScaling <- 100 # cm^2 per m^2

# Get time and location of individuals of interest
filtered_species_distr <- summary_shift %>%
  select(SpeciesPool, SpeciesID, last_year_alive, Shift, diff) %>%
  right_join(species_distr, ., by = c("SpeciesPool", "SpeciesID")) %>%
  filter(Year <= last_year_alive & Year >= last_year_alive - 20) %>%
  mutate(AreaOccupied = (Mass^(2/3)) / SurfaceBiomassScaling) # Area needed for individual

# Get species traits including FN for potential voxel counting below
species_traits <- filtered_species_distr %>%
  select(SpeciesID, SpeciesPool, AreaOccupied) %>%
  group_by(SpeciesID, SpeciesPool) %>%
  summarize(MinAreaOccupied = min(AreaOccupied, na.rm = TRUE),
            MaxAreaOccupied = max(AreaOccupied, na.rm = TRUE),
            MeanAreaOccupied = mean(AreaOccupied, na.rm = TRUE)
  ) %>%
  ungroup(.) %>%
  left_join(niches, ., by = c("SpeciesPool", "SpeciesID")) %>%
  mutate(MinLightScaled = MinLight / 900,
         MaxLightScaled = MaxLight / 900)

# - 2. Loop over years and extract MC matrices and env. conditions
scenario_dir <- list("CC" = "climdata_era5_cmip6_1981-2100_ssp245",
                     "No CC" = "climdata_era5_cmip6_1981-2100_ssp245_no_cc")
species_env <- NA
RealizedPotentialSpace <- NA
for (scenario in c("CC", "No CC")) {
  print(scenario)
  for (forest in 0:2) {
    print(forest)
    species_forest <- filtered_species_distr %>% filter(ForestID == forest, Scenario == scenario)

    for (ts in sort(unique(species_forest$TimeStep))) {
      print(ts)
      print("------")
      species_locations <- species_forest %>% filter(TimeStep == ts)

      # Load mC matrix
      mh <- readRDS(file.path(base_dir, scenario_dir[[scenario]], "a1_2",
                              paste0("forest", forest), paste0("microhabitatMatrix", ts, ".rds")))
      mc <- mh[,,, c(3,4,5)] # Get microclimate
      substrate_area <- mh[,,, 1] # Get substrate area

      # For each individual, get the max non-zero value in substrate_area[X, Y, ]
      species_locations$VegetationHeight <- mapply(function(x, y, areaNeeded) {
          layer_values <- substrate_area[x, y, ]
          max_non_zero <- max(which(layer_values >= areaNeeded))
        }, species_locations$X, species_locations$Y, species_locations$AreaOccupied)

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

      # -- Get the number of occupied Voxels per species
      occupiedSpace <- species_locations %>%
        group_by(SpeciesPool, SpeciesID) %>%
        summarize(nOccupiedVoxels = n_distinct(X, Y, Z),
                  nOccupiedHeights = n_distinct(Z),
                  nOccupiedColumns = n_distinct(X, Y),
                  .groups = "drop")

      # -- Get the number of suitable voxels per species
      currSpeciesTraits <- species_traits %>%
        semi_join(species_locations, by = c("SpeciesID", "SpeciesPool")) %>%
        distinct(SpeciesID, SpeciesPool, .keep_all = TRUE)

      for (i in seq_len(nrow(currSpeciesTraits))) {
        # Compute number of suitable heights/voxels
        result <- compute_suitable_spaces(
          mc, substrate_area, currSpeciesTraits, i
        )

        currSpeciesTraits[i, "nSuitableHeights"] <- result$nSuitableHeights
        currSpeciesTraits[i, "nSuitableVoxels"] <- result$nSuitableVoxels
        currSpeciesTraits[i, "nSuitableLight"] <- result$nSuitableLight
        currSpeciesTraits[i, "nSuitableHum"] <- result$nSuitableHum
        currSpeciesTraits[i, "nSuitableTemp"] <- result$nSuitableTemp
        currSpeciesTraits[i, "nSuitableArea"] <- result$nSuitableArea
        currSpeciesTraits[i, "nSuitableColumns"] <- result$nSuitableColumns
      }

      # Merge the two and add all meta info
      currRealizedPotentialSpace <- currSpeciesTraits %>%
        select(nSuitableColumns, nSuitableHeights, nSuitableVoxels, nSuitableLight, nSuitableHum, nSuitableTemp,
               nSuitableArea, SpeciesID, SpeciesPool) %>%
        left_join(occupiedSpace, ., by = c("SpeciesID", "SpeciesPool")) %>%
        mutate(TimeStep = ts, Forest = forest, Scenario = scenario)

      if (all(is.na(species_env))) {
        species_env <- species_locations
        RealizedPotentialSpace <- currRealizedPotentialSpace
      } else {
        species_env <- rbind(species_locations, species_env)
        RealizedPotentialSpace <- rbind(currRealizedPotentialSpace, RealizedPotentialSpace)
      }
    }
  }
}

indiv_env_niches <- left_join(species_env, niches, by = c("SpeciesID", "SpeciesPool"))
write_csv(indiv_env_niches, file.path(base_dir, "individuals_cc_vs_no_cc_niches_env_cond_shift.csv"))

# -- Summarize suitable/occupied voxels/heights as potential geographical range filling

GeoPotRangeFilling <- RealizedPotentialSpace %>%
  group_by(Scenario, SpeciesPool, SpeciesID) %>%
  summarize(
    GeoPotAreaRangeFill = (sum(nOccupiedVoxels, na.rm = TRUE) / sum(nSuitableArea, na.rm = TRUE)) * 100,
    GeoPotTempRangeFill = (sum(nOccupiedVoxels, na.rm = TRUE) / sum(nSuitableTemp, na.rm = TRUE)) * 100,
    GeoPotHumRangeFill = (sum(nOccupiedVoxels, na.rm = TRUE) / sum(nSuitableHum, na.rm = TRUE)) * 100,
    GeoPotLightRangeFill = (sum(nOccupiedVoxels, na.rm = TRUE) / sum(nSuitableLight, na.rm = TRUE)) * 100,
    GeoPotRangeFill = (sum(nOccupiedVoxels, na.rm = TRUE) / sum(nSuitableVoxels, na.rm = TRUE)) * 100,
    VerticalPotRangeFill = (sum(nOccupiedHeights, na.rm = TRUE) / sum(nSuitableHeights, na.rm = TRUE)) * 100,
    Geo2DRangeFill  = (sum(nOccupiedColumns) / sum(nSuitableColumns)) * 100,
  ) %>%
  ungroup(.)

# -- Summarize env. conditions and potential conditions

species_env_agg <- species_env %>%
  group_by(Scenario, SpeciesPool, SpeciesID, Shift) %>%
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

# Merge realized and fundamental niches and geographical range filling
species_env_niches <- left_join(species_env_agg, niches, by = c("SpeciesID", "SpeciesPool"))

# -- Constrain realizable values by potential niche

species_env_niches_c <- species_env_niches %>% # Constrain potential limits with fundamental niche
  mutate(
    MaxPotTemp = pmin(MaxPotTemp, MaxTemp),
    MinPotTemp = pmax(MinPotTemp, MinTemp),
    MaxPotHum = pmin(MaxPotHum, MaxHum),
    MinPotHum = pmax(MinPotHum, MinHum),
    MaxPotLight = pmin(MaxPotLight, MaxLight),
    MinPotLight = pmax(MinPotLight, MinLight),
  )

# -- Compute niche filling and range size

range_filling <- species_env_niches_c %>%
  mutate(TempFNFill = ((MaxTempCond - MinTempCond) / (MaxTemp - MinTemp)) * 100,
         HumFNFill = ((MaxHumCond - MinHumCond) / (MaxHum - MinHum)) * 100,
         LightFNFill = ((MaxLightCond - MinLightCond) / (MaxLight - MinLight)) * 100,
         #CombFNFill = (TempFNFill * HumFNFill * LightFNFill) / 10000,
         TempPNSize = MaxPotTemp - MinPotTemp,
         HumPNSize = MaxPotHum - MinPotHum,
         LightPNSize = MaxPotLight - MinPotLight,
         TempPNFill = ((MaxTempCond - MinTempCond) / TempPNSize) * 100,
         HumPNFill = ((MaxHumCond - MinHumCond) / HumPNSize) * 100,
         LightPNFill = ((MaxLightCond - MinLightCond) / LightPNSize) * 100,
  )

# Merge with geographical range filling table
niche_range_filling <- inner_join(range_filling, GeoPotRangeFilling,
                                            by = c("Scenario", "SpeciesPool", "SpeciesID"))

# Reshape data to long format for ggplot
range_filling_long <- niche_range_filling %>%
  select(Scenario, TempFNFill, HumFNFill, LightFNFill, TempPNFill, HumPNFill, LightPNFill, GeoPotRangeFill,
         Geo2DRangeFill, VerticalPotRangeFill, ends_with("Size")) %>%
  pivot_longer(
    cols = c(TempFNFill, HumFNFill, LightFNFill, TempPNFill, HumPNFill, LightPNFill, GeoPotRangeFill, Geo2DRangeFill,
             VerticalPotRangeFill, ends_with("Size")),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      TempFNFill = "Temperature FN filling",
                      HumFNFill = "Relative humidity FN filling",
                      LightFNFill = "Light FN filling",
                      TempPNFill = "Temperature PN filling",
                      HumPNFill = "Relative humidity PN filling",
                      LightPNFill = "Light PN filling",
                      TempPNSize = "Temperature PN size",
                      HumPNSize = "Relative humidity PN size",
                      LightPNSize = "Light PN size",
                      GeoPotRangeFill = "3D range filling",
                      Geo2DRangeFill = "Horizontal range filling",
                      VerticalPotRangeFill = "Vertical range filling"
    )
  )

# Add alternative scenario names
range_filling_long$ScenarioNames <- factor(
  range_filling_long$Scenario,
  levels = c("No CC", "CC"),
  labels = c("Baseline", "Climate change")
)

write_csv(range_filling_long, file.path(base_dir, "niche_range_filling_cc_vs_no_cc.csv"))

# ----------------------------------- Relating shift and niche filling differences ----------------------------------- #

# - Niche filling difference vs. shift

range_diff <- range_filling %>%
  select(SpeciesID, SpeciesPool, Scenario, AvgDiff, TempFNFill, TempPNFill,
         HumFNFill, HumPNFill, LightFNFill, LightPNFill, MaxTemp, MaxHum, MinLight) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = c(TempFNFill, TempPNFill, HumFNFill, HumPNFill, LightFNFill,
                    LightPNFill, MaxTemp, MaxHum, MinLight),
    names_sep = "_"
  ) %>%
  mutate(
    # Compute differences (CC - NoCC)
    Diff_TempFNFill = TempFNFill_CC - `TempFNFill_No CC`,
    Diff_TempPNFill = TempPNFill_CC - `TempPNFill_No CC`,
    Diff_HumFNFill = HumFNFill_CC - `HumFNFill_No CC`,
    Diff_HumPNFill = HumPNFill_CC - `HumPNFill_No CC`,
    Diff_LightFNFill = LightFNFill_CC - `LightFNFill_No CC`,
    Diff_LightPNFill = LightPNFill_CC - `LightPNFill_No CC`,
  )

# Compute correlartions
cor.test(range_diff$Diff_TempFNFill, range_diff$AvgDiff, use = "complete.obs")
cor.test(range_diff$Diff_HumFNFill, range_diff$AvgDiff, use = "complete.obs")
cor.test(range_diff$Diff_LightFNFill, range_diff$AvgDiff, use = "complete.obs")
cor.test(range_diff$MaxTemp_CC, range_diff$AvgDiff, use = "complete.obs")

# ---------- Scatter plot of shift vs niche filling difference ---------- #

# - Temperature
pdf(file.path(DirectoryPlots, "scatter_shift_vs_niche_fill_diff_temp.pdf"), width = 7, height = 6)
ggplot(range_diff, aes(x = Diff_TempFNFill, y = AvgDiff)) +
  geom_point(size = 2, alpha = 0.7, color = "#004aad") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_smooth(method = "lm", se = TRUE, color = "#f7766e") +
  labs(
    x = expression(Delta~"Temperature niche filling (%)"~(CC - No~CC)),
    y = expression(Delta~"Position shift (m)"~(CC - No~CC))
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()

# - Humidity
pdf(file.path(DirectoryPlots, "scatter_shift_vs_niche_fill_diff_hum.pdf"), width = 7, height = 6)
ggplot(range_diff, aes(x = Diff_HumFNFill, y = AvgDiff)) +
  geom_point(size = 2, alpha = 0.7, color = "#004aad") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_smooth(method = "lm", se = TRUE, color = "#f7766e") +
  labs(
    x = expression(Delta~"Relatve humidity niche filling (%)"~(CC - No~CC)),
    y = expression(Delta~"Position shift (m)"~(CC - No~CC))
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()

# - Light
pdf(file.path(DirectoryPlots, "scatter_shift_vs_niche_fill_diff_light.pdf"), width = 7, height = 6)
ggplot(range_diff, aes(x = Diff_LightFNFill, y = AvgDiff)) +
  geom_point(size = 2, alpha = 0.7, color = "#004aad") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_smooth(method = "lm", se = TRUE, color = "#f7766e") +
  labs(
    x = expression(Delta~"Light niche filling (%)"~(CC - No~CC)),
    y = expression(Delta~"Position shift (m)"~(CC - No~CC))
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()


# ---------------------------------------- Compare range sizes ---------------------------------------- #

# --- Compute range sizes

range_sizes <- species_env_niches_c %>%
  mutate(TempRealizedSize = MaxTempCond - MinTempCond,
         HumRealizedSize = MaxHumCond - MinHumCond,
         LightRealizedSize = MaxLightCond - MinLightCond,
         TempRangeSize = MaxPotTemp - MinPotTemp,
         HumRangeSize = MaxPotHum - MinPotHum,
         LightRangeSize = MaxPotLight - MinPotLight,
         TempNicheSize = MaxTemp - MinTemp,
         HumNicheSize = MaxHum - MinHum,
         LightNicheSize = MaxLight - MinLight) %>%
  select(Scenario, SpeciesPool, SpeciesID, TempRealizedSize, HumRealizedSize, LightRealizedSize, TempRangeSize, HumRangeSize, LightRangeSize,
         TempNicheSize, HumNicheSize, LightNicheSize) %>%
  pivot_longer(
    cols = c(TempRealizedSize, HumRealizedSize, LightRealizedSize, TempRangeSize, HumRangeSize, LightRangeSize,
             TempNicheSize, HumNicheSize, LightNicheSize),
    names_to = "Variable",
    values_to = "Size"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      TempRealizedSize = "Temperature realized range",
                      HumRealizedSize = "Relative humidity realized range",
                      LightRealizedSize = "Light realized range",
                      TempRangeSize = "Temperature potential range",
                      HumRangeSize = "Relative humidity potential range",
                      LightRangeSize = "Light potential range",
                      TempNicheSize = "Temperature fundamental niche",
                      HumNicheSize = "Relative humidity fundamental niche",
                      LightNicheSize = "Light fundamental niche"
    )
  )

# - Plot

# Create boxplots
pLight <- ggplot(range_sizes %>% filter(grepl('Light', Variable)),
                 aes(x = Variable, y = Size, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  labs(
    x = "",
    y = "Light range size",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
pTemp <- ggplot(range_sizes %>% filter(grepl('Temp', Variable)),
                 aes(x = Variable, y = Size, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  labs(
    x = "",
    y = "Temperature range size",
    fill = "Scenario"
  )
pHum <- ggplot(range_sizes %>% filter(grepl('Relative', Variable)),
                 aes(x = Variable, y = Size, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  labs(
    x = "",
    y = "Relative humidity range size",
    fill = "Scenario"
  )

output_file3 <- file.path(DirectoryPlots, "boxplot_range_niche_realized_sizes_v2.pdf")
pdf(output_file3, width = 7, height = 6)
combined_plot <-
  (pTemp + pHum + pLight) +
  plot_layout(ncol = 1) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1.5, "cm")
  )
print(combined_plot)
dev.off()