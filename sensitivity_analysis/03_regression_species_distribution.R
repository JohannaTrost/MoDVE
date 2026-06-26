library(dplyr)
library(readr)
library(lmerTest)
library(lattice)
library(car)
library(performance)

DirectoryModelResults <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")

fileName <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_SpeciesVertical_mcGradients.csv"

verticalDistr <- read_csv(file.path(DirectoryModelResults, fileName), show_col_types = FALSE)

# - Regression analysis

# 1. Perform mixed-effects regression for mean position

verticalDistr <- verticalDistr %>%
  mutate(
    timeStep_z = as.numeric(scale(timeStep)),
  )
mem_mean_hgt <- lmer(meanHeight ~ scenario * timeStep + (1|speciesPool),
                          data = verticalDistr)
summary(mem_mean_hgt)
r2(mem_mean_hgt)

# Without time component

lm_mean_hgt_ts <- verticalDistr %>%
  filter(timeStep == max(verticalDistr$timeStep)) %>%
  lm(meanHeight ~ scenario, data = .)

summary(lm_mean_hgt_ts)

# 2. For breadth of distribution (Variance)

mem_mean_var <- lmer(varHeight ~ scenario * timeStep + (1|speciesPool),
                          data = verticalDistr)
summary(mem_mean_var)
r2(mem_mean_var)

# Without time component

lm_var_hgt_ts <- verticalDistr %>%
  filter(timeStep == max(verticalDistr$timeStep)) %>%
  lm(varHeight ~ scenario, data = .)

summary(lm_var_hgt_ts)

# 3. For breadth of distribution (IQR)

mem_mean_iqr <- lmer(iqrHeight ~ scenario * timeStep + (1|speciesPool),
                     data = verticalDistr)
summary(mem_mean_iqr)
r2(mem_mean_iqr)
# Without time component
lm_iqr_hgt_ts <- verticalDistr %>%
  filter(timeStep == max(verticalDistr$timeStep)) %>%
  lm(iqrHeight ~ scenario, data = .)
summary(lm_iqr_hgt_ts)