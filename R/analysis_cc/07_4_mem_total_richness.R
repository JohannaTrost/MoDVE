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
richness <- read_csv(file.path(base_dir, "a5_diversity_cc_vs_no_cc.csv"))

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc")

# Sort scenarios and arrange data
richness_filtered <- richness %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario) %>%
  filter(Year >= 2020) # Remove years of burnin

# z-score year
richness_filtered$Year_c <- (richness_filtered$Year - mean(richness_filtered$Year)) / sd(richness_filtered$Year)

# Check distribution of responses
pdf(file.path(DirectoryPlots, "total_richness_hist.pdf"))
hist(richness_filtered$Richness)
dev.off() # -> Poisson

pdf(file.path(DirectoryPlots, "total_abundance_hist.pdf"))
hist(richness_filtered$Abundance)
dev.off() # -> Poisson