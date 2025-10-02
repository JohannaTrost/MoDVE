library(readr)
library(lme4)
library(nlme)
library(lmerTest)
library(dplyr)
library(tidyr)

# Repsonse vars:
# - Avg. species position
# - Species range size
# - Species richness (peak)
# - Species abundance (peak)
# - Species diversity (peak)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")
species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

# Get stats
species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            .groups = "drop")

# Compute CC position shift between scenarios at species level for average of last 20 years
# (neg. = donwshift with CC)
diff_stats_species_level <- species_distr_stats %>%
  filter(Year >= 2079) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = c(Position, Mass)
  ) %>%
  mutate(
    Mass_diff = Mass_CC - `Mass_No CC`,
    Position_diff = Position_CC - `Position_No CC`
  ) %>%
  dplyr::select(SpeciesPool, SpeciesID, ForestID, Year, Mass_diff, Position_diff) %>%
  group_by(SpeciesPool, SpeciesID, ForestID) %>%
  summarise(
    Mass_diff = mean(Mass_diff, na.rm = TRUE),
    Position_diff = mean(Position_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
      Mass_diff = mean(Mass_diff, na.rm = TRUE),
      Mass_diff_sd = sd(Mass_diff, na.rm = TRUE),
      Position_diff = mean(Position_diff, na.rm = TRUE),
      Position_diff_sd = sd(Position_diff, na.rm = TRUE),
      .groups = "drop"
  )

diff_stats_species_level %>%
  filter(!is.na(Position_diff)) %>%
  write_csv(., file.path(base_dir, "a5_species_level_position_mass_shift_2080-2100.csv"))

# Compute CC position shift between scenarios aggregated at species pool level for average of last 20 years
diff_stats_species_level %>%
  group_by(SpeciesPool) %>%
    summarize(
        mean_pos_diff = mean(Position_diff, na.rm = TRUE),
        sd_pos_diff = sd(Position_diff, na.rm = TRUE),
        mean_mass_diff = mean(Mass_diff, na.rm = TRUE),
        sd_mass_diff = sd(Mass_diff, na.rm = TRUE),
        .groups = "drop"
    ) %>%
  bind_rows(
    tibble(
      SpeciesPool = NA,
      mean_pos_diff = mean(.$mean_pos_diff),
      sd_pos_diff = sd(.$mean_pos_diff),
      mean_mass_diff = mean(.$mean_mass_diff),
      sd_mass_diff = sd(.$mean_mass_diff)
    )
  ) %>%
  write_csv(., file.path(base_dir, "a5_species_pool_level_position_mass_shift_2080-2100.csv"))