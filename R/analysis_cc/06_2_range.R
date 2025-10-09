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
  summarize(IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop")

# Compute CC range shift between scenarios at species level for average of last 20 years
# (neg. = donwshift with CC)
diff_stats_species_level <- species_distr_stats %>%
  filter(Year >= 2079) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = c(IQR, Range)
  ) %>%
  mutate(
    IQR_diff = IQR_CC - `IQR_No CC`,
    Range_diff = Range_CC - `Range_No CC`
  ) %>%
  dplyr::select(SpeciesPool, SpeciesID, ForestID, Year, IQR_diff, Range_diff) %>%
  group_by(SpeciesPool, SpeciesID, ForestID) %>%
  summarise(
    IQR_diff = mean(IQR_diff, na.rm = TRUE),
    Range_diff = mean(Range_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
      IQR_diff = mean(IQR_diff, na.rm = TRUE),
      IQR_diff_sd = sd(IQR_diff, na.rm = TRUE),
      Range_diff = mean(Range_diff, na.rm = TRUE),
      Range_diff_sd = sd(Range_diff, na.rm = TRUE),
      .groups = "drop"
  )

diff_stats_species_level %>%
  filter(!is.na(IQR_diff)) %>%
  write_csv(., file.path(base_dir, "a5_species_level_range_iqr_shift_2080-2100.csv"))

# Compute CC position shift between scenarios aggregated at species pool level for average of last 20 years
diff_stats_species_level %>%
  group_by(SpeciesPool) %>%
    summarize(
        mean_Range_diff = mean(Range_diff, na.rm = TRUE),
        sd_Range_diff = sd(Range_diff, na.rm = TRUE),
        mean_IQR_diff = mean(IQR_diff, na.rm = TRUE),
        sd_IQR_diff = sd(IQR_diff, na.rm = TRUE),
        .groups = "drop"
    ) %>%
  bind_rows(
    tibble(
      SpeciesPool = NA,
      mean_Range_diff = mean(.$mean_Range_diff),
      sd_Range_diff = sd(.$mean_Range_diff),
      mean_IQR_diff = mean(.$mean_IQR_diff),
      sd_IQR_diff = sd(.$mean_IQR_diff)
    )
  ) %>%
  write_csv(., file.path(base_dir, "a5_species_pool_level_range_iqr_shift_2080-2100.csv"))