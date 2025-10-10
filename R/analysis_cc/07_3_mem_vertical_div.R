library(readr)
library(lme4)
library(nlme)
library(lmerTest)
library(dplyr)
library(merTools)
library(DHARMa)
library(lattice)
library(ggplot2)
library(patchwork)  # for combining plots
library(glmmTMB)
library(lattice)
library(DHARMa)
library(MASS)
library(performance)
library(tidyr)
library(purrr)

# Repsonse vars:
# - Avg. species position
# - Species range size
# - Species richness (peak)
# - Species abundance (peak)
# - Species diversity (peak)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")
div_peak <- read_csv(file.path(base_dir, "a5_vertical_diversity_peak_cc_vs_no_cc.csv"))

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/VerticalDiversity")
dir.create(DirectoryPlots)

# Sort scenarios and arrange data
div_peak <- div_peak %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

div_peak_filtered <- div_peak %>%
  filter(Year >= 2020) %>%
  # Identify all Scenario-Year-ForestID-SpeciesPool combos where Richness < 2
  anti_join(
    div_peak %>%
      filter(level_5 == "Richness", Count < 2) %>%
      dplyr::select(Scenario, Year, ForestID, SpeciesPool),
    by = c("Scenario", "Year", "ForestID", "SpeciesPool")
  )

div_peak_filtered$Year_c <- (div_peak_filtered$Year - mean(div_peak_filtered$Year)) / sd(div_peak_filtered$Year)

