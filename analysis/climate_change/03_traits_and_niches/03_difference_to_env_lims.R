# -----
# Compute distances of species' environmental conditions with max. tolerance limits and optima
# (for temperature, relative humidity and light)

library(readr)
library(dplyr)

base_dir <- file.path("../modve_data/modve_output/regua")
indiv_env_niches <- read_csv(file.path(base_dir, "individuals_cc_vs_no_cc_niches_env_cond_shift.csv"))

# - Compute distance to tolerances
distances <- indiv_env_niches %>%
  mutate(MaxTempDist = abs(MaxTemp - Temp),
         MaxHumDist = abs(MaxHum - Hum),
         MinLightDist = abs(MinLight - Light),
         OptTempDist = abs(OptimumTemp - Temp),
         OptHumDist = abs(OptimumHum - Hum),
         OptLightDist = abs(OptimumLight - Light)
  ) %>%
  group_by(Scenario, SpeciesPool, SpeciesID) %>%
  summarise(MaxTempDist = mean(MaxTempDist, na.rm = TRUE),
            MaxHumDist = mean(MaxHumDist, na.rm = TRUE),
            MinLightDist = mean(MinLightDist, na.rm = TRUE),
            OptTempDist = mean(OptTempDist, na.rm = TRUE),
            OptHumDist = mean(OptHumDist, na.rm = TRUE),
            OptLightDist = mean(OptLightDist, na.rm = TRUE),
  )

# - Distance to max. limit

# Reshape data to long format for ggplot
max_dists_long <- distances %>%
  select(Scenario, SpeciesPool, SpeciesID, MaxTempDist, MaxHumDist, MinLightDist) %>%
  pivot_longer(
    cols = c(MaxTempDist, MaxHumDist, MinLightDist),
    names_to = "Variable",
    values_to = "DistanceToMax"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      MaxTempDist = "Temperature",
                      MaxHumDist = "Relative humidity",
                      MinLightDist = "Light")
  )

# - Distance to optimum max. limit

# Reshape data to long format for ggplot
opt_dists_long <- distances %>%
  select(Scenario, SpeciesPool, SpeciesID, OptTempDist, OptHumDist, OptLightDist) %>%
  pivot_longer(
    cols = c(OptTempDist, OptHumDist, OptLightDist),
    names_to = "Variable",
    values_to = "DistanceToOpt"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      OptTempDist = "Temperature",
                      OptHumDist = "Relative humidity",
                      OptLightDist = "Light")
  )

# Merge and save distance df
dist_df <- inner_join(opt_dists_long, max_dists_long, by = c("Scenario", "SpeciesPool", "SpeciesID", "Variable"))
write_csv(dist_df, file.path(base_dir, "opt_max_distances_cc_vs_no_cc.csv"))