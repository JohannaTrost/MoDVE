library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)
library(vegan)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

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

year <- 2080
species_shift <- species_distr_stats %>%
  filter(Year >= year) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    AvgDiff = mean(diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Shift = ifelse(AvgDiff > 0, "Upward", "Downward")) %>%
  drop_na(.)