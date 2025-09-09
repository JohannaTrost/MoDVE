library(dplyr)
library(readr)
library(tidyverse)

###############################################################################################
#                                        Diversity data                                       #
###############################################################################################

DirectoryModelResults <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")

divFile1 <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxRichness_mcGradients.csv"
divFile2 <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxShannonDiv_mcGradients.csv"

maxRichness <- read_csv(file.path(DirectoryModelResults, divFile1), show_col_types = FALSE)
maxShannon <- read_csv(file.path(DirectoryModelResults, divFile2), show_col_types = FALSE)

# - Compute relative changes compared to gradient scenarios

############# Richness peak position #############

relative_change_richness <- maxRichness %>%
  # Reshape so each scenario has its own column for maxZ
  select(X, Y, speciesPool, timeStep, scenario, maxZ) %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = maxZ,
    names_prefix = "scenario_"
  ) %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = (scenario_0   - scenario_1) / scenario_1,
    relChange_0.5 = (scenario_0.5 - scenario_1) / scenario_1,
    relChange_1.5 = (scenario_1.5 - scenario_1) / scenario_1
  )

# - Summarize relative changes across space and time

relative_change_richness_summary <- relative_change_richness %>%
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
relative_change_richness_summary

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
      ~sd(.x, na.rm = TRUE),
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

############# Richness peak value #############

relative_change_max_richness <- maxRichness %>%
  rename(maxR = maxRichness) %>%
  # Reshape so each scenario has its own column for maxZ
  select(X, Y, speciesPool, timeStep, scenario, maxR) %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = maxR,
    names_prefix = "scenario_"
  ) %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = (scenario_0   - scenario_1) / scenario_1,
    relChange_0.5 = (scenario_0.5 - scenario_1) / scenario_1,
    relChange_1.5 = (scenario_1.5 - scenario_1) / scenario_1
  )

# - Summarize relative changes across space and time

relative_change_max_richness_summary <- relative_change_max_richness %>%
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
relative_change_max_richness_summary

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
      ~sd(.x, na.rm = TRUE),
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

############### Shannon diversity peak position #############

relative_change_div <- maxShannon %>%
  # Reshape so each scenario has its own column for maxZ
  select(X, Y, speciesPool, timeStep, scenario, maxZ) %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = maxZ,
    names_prefix = "scenario_"
  ) %>%
  # Compute relative change compared to baseline scenario==1
  mutate(
    relChange_0   = (scenario_0   - scenario_1) / scenario_1,
    relChange_0.5 = (scenario_0.5 - scenario_1) / scenario_1,
    relChange_1.5 = (scenario_1.5 - scenario_1) / scenario_1
  )

# - Summarize relative changes across space and time

relative_change_div_summary <- relative_change_div %>%
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
relative_change_div_summary

# Compute mean and sd of relative changes across species pools
relative_change_div_summary %>%
  summarise(
    across(
      starts_with("meanRelChange_"),
      ~mean(.x, na.rm = TRUE),
      .names = "overallMean_{.col}"
    ),
    across(
      starts_with("meanRelChange_"),
      ~sd(.x, na.rm = TRUE),
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

###############################################################################################
#                             Species distribution data                                       #
###############################################################################################

DirectoryModelResults <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")

fileName <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_SpeciesVertical_mcGradients.csv"

verticalDistr <- read_csv(file.path(DirectoryModelResults, fileName), show_col_types = FALSE)

############### Species Position #############

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
    relChange_0   = (scenario_0   - scenario_1) / scenario_1,
    relChange_0.5 = (scenario_0.5 - scenario_1) / scenario_1,
    relChange_1.5 = (scenario_1.5 - scenario_1) / scenario_1
  )

# - Summarize relative changes across space and time

relative_change_distr_summary <- relative_change_distr %>%
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
relative_change_distr_summary

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
      ~sd(.x, na.rm = TRUE),
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

############### Species range (Variance) #############

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
relative_change_range_summary

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


############### Species range (IQR) #############

relative_change_iqr <- verticalDistr %>%
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
relative_change_iqr_summary

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