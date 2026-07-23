# -----
# Analyse relative changes in species position, range and richness across MC gradinet steepness scenarios

library(dplyr)
library(readr)
library(tidyverse)

base_dir <- file.path("../modve_data/modve_output/pirineus/scenarios")
maxRichness <- read_csv(file.path(base_dir, "vertical_richness_steepness.csv"))

if (!file.exists(file.path(base_dir, "vertical_richness_steepness.csv"))) {
  stop(file.path(base_dir, "vertical_richness_steepness.csv"), " missing. Compute vertical richness for sensitivity ",
       "analysis using the script analysis/sensitivity/02_compute_richness.py")
}

# ------------------------------------ Richness peak position ------------------------------------ #

relative_change_richness <- maxRichness %>%
  # Reshape so each scenario has its own column for maxZ
  select(SpeciesPool, Year, Scenario, maxHeight) %>%
  tidyr::pivot_wider(
    names_from = Scenario,
    values_from = maxHeight,
    names_prefix = "scenario_"
  ) %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = 100*((scenario_0   - scenario_1) / scenario_1),
    relChange_0.5 = 100*((scenario_0.5 - scenario_1) / scenario_1),
    relChange_1.5 = 100*((scenario_1.5 - scenario_1) / scenario_1)
  )

# - Summarize relative changes across space and time

relative_change_richness_summary <- relative_change_richness %>%
  group_by(SpeciesPool) %>%
    summarise(
        meanRelChange_0   = mean(relChange_0, na.rm = TRUE),
        meanRelChange_0.5 = mean(relChange_0.5, na.rm = TRUE),
        meanRelChange_1.5 = mean(relChange_1.5, na.rm = TRUE),
        sdRelChange_0     = sd(relChange_0, na.rm = TRUE),
        sdRelChange_0.5   = sd(relChange_0.5, na.rm = TRUE),
        sdRelChange_1.5   = sd(relChange_1.5, na.rm = TRUE),
        nVoxels           = n()
    )
print(relative_change_richness_summary)

# Compute mean and sd of relative changes across species pools
relative_change_richness_summary %>%
  summarise(
    across(
      starts_with("meanRelChange_"),
      ~mean(.x, na.rm = TRUE),
      .names = "overallMean_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.95, na.rm = TRUE),
      .names = "overallp95_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.05, na.rm = TRUE),
      .names = "overallp5_{.col}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "Scenario"),
    names_pattern = "(overallMean|overallp95|overallp5)_meanRelChange_(.*)"
  ) %>%
  rename(
    overallMeanRelChange = overallMean,
    overallp95RelChange   = overallp95,
    overallp5RelChange   = overallp5
  )

# ------------------------------------ Richness peak value ------------------------------------ #

relative_change_max_richness <- maxRichness %>%
  rename(maxR = maxRichness) %>%
  # Reshape so each scenario has its own column for maxZ
  select(SpeciesPool, Year, Scenario, maxR) %>%
  tidyr::pivot_wider(
    names_from = Scenario,
    values_from = maxR,
    names_prefix = "scenario_"
  ) %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = 100*((scenario_0   - scenario_1) / scenario_1),
    relChange_0.5 = 100*((scenario_0.5 - scenario_1) / scenario_1),
    relChange_1.5 = 100*((scenario_1.5 - scenario_1) / scenario_1)
  )

# - Summarize relative changes across space and time

relative_change_max_richness_summary <- relative_change_max_richness %>%
  group_by(SpeciesPool) %>%
    summarise(
        meanRelChange_0   = mean(relChange_0, na.rm = TRUE),
        meanRelChange_0.5 = mean(relChange_0.5, na.rm = TRUE),
        meanRelChange_1.5 = mean(relChange_1.5, na.rm = TRUE),
        sdRelChange_0     = sd(relChange_0, na.rm = TRUE),
        sdRelChange_0.5   = sd(relChange_0.5, na.rm = TRUE),
        sdRelChange_1.5   = sd(relChange_1.5, na.rm = TRUE),
        nVoxels           = n()
    )
print(relative_change_max_richness_summary)

# Compute mean and sd of relative changes across species pools
relative_change_max_richness_summary %>%
  summarise(
    across(
      starts_with("meanRelChange_"),
      ~mean(.x, na.rm = TRUE),
      .names = "overallMean_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.95, na.rm = TRUE),
      .names = "overallp95_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.05, na.rm = TRUE),
      .names = "overallp5_{.col}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "Scenario"),
    names_pattern = "(overallMean|overallp95|overallp5)_meanRelChange_(.*)"
  ) %>%
  rename(
    overallMeanRelChange = overallMean,
    overallp5RelChange   = overallp5,
    overallp95RelChange   = overallp95
  )

# ------------------------------------ Species Position ------------------------------------ #

fileName <- "SpeciesVertical_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_mcGradients.csv"
verticalDistr <- read_csv(file.path(base_dir, fileName), show_col_types = FALSE)

relative_change_distr <- verticalDistr %>%
  # Reshape so each scenario has its own column for meanHeight
  select(SpeciesID, speciesPool, timeStep, scenario, meanHeight) %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = meanHeight,
    names_prefix = "scenario_"
  ) %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = 100*((scenario_0   - scenario_1) / scenario_1),
    relChange_0.5 = 100*((scenario_0.5 - scenario_1) / scenario_1),
    relChange_1.5 = 100*((scenario_1.5 - scenario_1) / scenario_1)
  )

# Check if normal
shapiro.test(sample(relative_change_distr$scenario_0, 1000))
shapiro.test(sample(relative_change_distr$scenario_0.5, 1000))
shapiro.test(sample(relative_change_distr$scenario_1.5, 1000))
shapiro.test(sample(relative_change_distr$scenario_1, 1000))

# - Summarize relative changes across space and time

relative_change_distr_summary <- relative_change_distr %>%
  filter(timeStep >= 187 & timeStep <= 197) %>%
  group_by(speciesPool) %>%
    summarise(
        meanRelChange_0   = mean(relChange_0, na.rm = TRUE),
        meanRelChange_0.5 = mean(relChange_0.5, na.rm = TRUE),
        meanRelChange_1.5 = mean(relChange_1.5, na.rm = TRUE),
        sdRelChange_0     = sd(relChange_0, na.rm = TRUE)* 1.96,
        sdRelChange_0.5   = sd(relChange_0.5, na.rm = TRUE)* 1.96,
        sdRelChange_1.5   = sd(relChange_1.5, na.rm = TRUE)* 1.96,
        nVoxels           = n()
    )
print(relative_change_distr_summary)

# Compute mean and sd of relative changes across species pools
relative_change_distr_summary %>%
  summarise(
    across(
      starts_with("meanRelChange_"),
      ~mean(.x, na.rm = TRUE),
      .names = "overallMean_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.05, na.rm = TRUE),
      .names = "overallp5_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.95, na.rm = TRUE),
      .names = "overallp95_{.col}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "Scenario"),
    names_pattern = "(overallMean|overallp5|overallp95)_meanRelChange_(.*)"
  ) %>%
  rename(
    overallMeanRelChange = overallMean,
    overallp5RelChange   = overallp5,
    overallp95RelChange   = overallp95
  )

# ------------------------------------ Species range (Variance) ------------------------------------ #

relative_change_range <- verticalDistr %>%
  # Reshape so each scenario has its own column for meanHeight
  select(SpeciesID, speciesPool, timeStep, scenario, varHeight) %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = varHeight,
    names_prefix = "scenario_"
  ) %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = (scenario_0   - scenario_1) / scenario_1,
    relChange_0.5 = (scenario_0.5 - scenario_1) / scenario_1,
    relChange_1.5 = (scenario_1.5 - scenario_1) / scenario_1
  )

# - Summarize relative changes across space and time

relative_change_range_summary <- relative_change_range %>%
  filter(timeStep >= 187 & timeStep <= 197) %>%
  group_by(speciesPool) %>%
    summarise(
        meanRelChange_0   = mean(relChange_0, na.rm = TRUE),
        meanRelChange_0.5 = mean(relChange_0.5, na.rm = TRUE),
        meanRelChange_1.5 = mean(relChange_1.5, na.rm = TRUE),
        sdRelChange_0     = sd(relChange_0, na.rm = TRUE),
        sdRelChange_0.5   = sd(relChange_0.5, na.rm = TRUE),
        sdRelChange_1.5   = sd(relChange_1.5, na.rm = TRUE),
        nVoxels           = n()
    )
print(relative_change_range_summary)

# Compute mean and sd of relative changes across species pools
relative_change_range_summary %>%
  mutate(across(starts_with("meanRelChange_"), ~ ifelse(is.infinite(.x), NA, .x))) %>%
  summarise(
    across(
      starts_with("meanRelChange_"),
      ~ mean(.x, na.rm = TRUE),
      .names = "overallMean_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~ sd(.x, na.rm = TRUE),
      .names = "overallSd_{.col}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "Scenario"),
    names_pattern = "(overallMean|overallSd)_meanRelChange_(.*)"
  ) %>%
  rename(
    overallMeanRelChange = overallMean,
    overallSdRelChange   = overallSd
  )

# ------------------------------------ Species range (IQR) ------------------------------------ #

relative_change_iqr <- verticalDistr %>%
  # Reshape so each scenario has its own column for meanHeight
  select(SpeciesID, speciesPool, timeStep, scenario, iqrHeight) %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = iqrHeight,
    names_prefix = "scenario_"
  ) %>%
  drop_na() %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = 100*((scenario_0   - scenario_1) / scenario_1),
    relChange_0.5 = 100*((scenario_0.5 - scenario_1) / scenario_1),
    relChange_1.5 = 100*((scenario_1.5 - scenario_1) / scenario_1)
  )

shapiro.test(sample(relative_change_iqr$scenario_0, 1000))
shapiro.test(sample(relative_change_iqr$scenario_0.5, 1000))
shapiro.test(sample(relative_change_iqr$scenario_1.5, 1000))
shapiro.test(sample(relative_change_iqr$scenario_1, 1000))

# - Summarize relative changes across space and time

relative_change_iqr_summary <- relative_change_iqr %>%
  filter(timeStep >= 187 & timeStep <= 197) %>%
  group_by(speciesPool) %>%
    summarise(
        meanRelChange_0   = mean(relChange_0, na.rm = TRUE),
        meanRelChange_0.5 = mean(relChange_0.5, na.rm = TRUE),
        meanRelChange_1.5 = mean(relChange_1.5, na.rm = TRUE),
        sdRelChange_0     = sd(relChange_0, na.rm = TRUE),
        sdRelChange_0.5   = sd(relChange_0.5, na.rm = TRUE),
        sdRelChange_1.5   = sd(relChange_1.5, na.rm = TRUE),
        nVoxels           = n()
    )
print(relative_change_iqr_summary)

# Compute mean and sd of relative changes across species pools
relative_change_iqr_summary %>%
  mutate(across(starts_with("meanRelChange_"), ~ ifelse(is.infinite(.x), NA, .x))) %>%
  summarise(
    across(
      starts_with("meanRelChange_"),
      ~ mean(.x, na.rm = TRUE),
      .names = "overallMean_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.05, na.rm = TRUE),
      .names = "overallp5_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~quantile(.x, 0.95, na.rm = TRUE),
      .names = "overallp95_{.col}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "Scenario"),
    names_pattern = "(overallMean|overallp95|overallp5)_meanRelChange_(.*)"
  ) %>%
  rename(
    overallMeanRelChange = overallMean,
    overallp5RelChange   = overallp5,
    overallp95RelChange   = overallp95
  )