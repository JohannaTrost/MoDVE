# -----
# Compute species position shifts between CC and baseline scenario
# across last 20 years of the species lifetime for all species that survive 2050
# (Assumes execution from project root)

library(readr)
library(dplyr)
library(tidyr)

base_dir <- file.path("../modve_data/modve_output/regua")
spec_distr_path <- file.path(base_dir, "species_distribution_cc_vs_no_cc.csv")
if (!file.exists(spec_distr_path)) {
  stop(
    spec_distr_path,
    " does not exist. Run analysis/climate_change/01_community/03_species_distributions.py ",
    "to generate the distribution file."
  )
}

# Load species distribution file
species_distr <- read_csv(spec_distr_path)

# Get stats
species_distr_stats <- species_distr %>%
  group_by(Scenario, ForestID, SpeciesPool, SpeciesID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop")

# Sort scenarios and arrange data
species_distr_stats <- species_distr_stats %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

# Compute position shift
shift <- species_distr_stats %>%
  dplyr::select(ForestID, Scenario, SpeciesPool, SpeciesID, TimeStep, Year, Position) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = Position
  ) %>%
  mutate(
    diff = `CC` - `No CC`
  ) %>%
  filter(!is.na(diff))

# Average species shift across forests for the last 20 years of each species that survives past 2050
summary_shift <- shift %>%
  group_by(SpeciesPool, SpeciesID) %>%
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    initial_diff = mean(.data$diff[.data$Year == 2000], na.rm = TRUE),
    valid_diff = .data$diff[!is.na(.data$diff)],
    diff = mean(tail(valid_diff, 20), na.rm = TRUE),
    diff_sem = sd(tail(valid_diff, 20), na.rm = TRUE) / sqrt(length(tail(valid_diff, 20))),
    last_year_alive = max(.data$Year[!is.na(.data$diff)]),
    .groups = "drop"
  ) %>%
  filter(last_year_alive >= 2050) %>%
  select(-valid_diff) %>%
  distinct()

# Determine upward shifting species
shift_species_q75 <- summary_shift %>%
  filter(diff >= quantile(diff, 0.75)) # Threshold is positive

# Save species shift stats and only upward shifting species
write_csv(summary_shift, file.path(base_dir, "species_shift_cc_vs_no_cc_last20_yrs_past_2050.csv"))
write_csv(shift_species_q75,
          file.path(base_dir, "species_shift_upwards_cc_vs_no_cc.csv"))