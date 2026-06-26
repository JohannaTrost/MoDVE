library(dplyr)
library(readr)
library(lmerTest)
library(lattice)
library(car)
library(performance)

DirectoryModelResults <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")

divFile1 <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxRichness_mcGradients.csv"
divFile2 <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxShannonDiv_mcGradients.csv"

maxRichness <- read_csv(file.path(DirectoryModelResults, divFile1), show_col_types = FALSE)
maxShannon <- read_csv(file.path(DirectoryModelResults, divFile2), show_col_types = FALSE)

# - Regression analysis

# - Perform mixed-effects regression for maxRichness

# MEM
maxRichness <- maxRichness %>%
  mutate(
    timeStep_z = as.numeric(scale(timeStep)),
    maxZ_z     = as.numeric(scale(maxZ))
  )
mem_richness_test <- lmer(maxRichness ~ scenario * timeStep_z * maxZ_z +
                          (1|X) + (1|Y) + (1|speciesPool),
                          data = maxRichness)
summary(mem_richness_test)

# R2
r2(mem_richness_test)

# --------------------------------------------------------------------------------------------------

# Avg across x, y
avg_richness <- maxRichness %>%
  group_by(scenario, speciesPool, timeStep_z, maxZ_z) %>%
  summarise(meanMaxRichness = mean(maxRichness), .groups = "drop")

# Regression on avg richness without x,y
mem_avg_richness <- lmer(meanMaxRichness ~ scenario * timeStep_z * maxZ_z + (1|speciesPool),
                          data = avg_richness)
summary(mem_avg_richness)
# R2
r2(mem_avg_richness)

# --------------------------------------------------------------------------------------------------

# Pick one Time step
mem_avg_richness_ts <-avg_richness %>%
  filter(timeStep_z == max(avg_richness$timeStep_z)) %>%
  lmer(meanMaxRichness ~ scenario * maxZ_z + (1|speciesPool), data = .)
summary(mem_avg_richness_ts)
# R2
r2(mem_avg_richness_ts)

# --------------------------------------------------------------------------------------------------

# Exploratory plot
dir.create("../../figs/sensitivity_analysis")
pdf("../../figs/sensitivity_analysis/maxRichness_scenario_speciesPool.pdf")
xyplot(maxRichness ~ factor(scenario) | factor(speciesPool), data = maxRichness)
dev.off()

# - Perfrom LLR to check if random effect is meaningful
# Full model (with all random effects)
mem_richness <- lmer(maxRichness ~ scenario * timeStep_z * maxZ_z +
                     (1|X) + (1|Y) + (1|speciesPool),
                     data = maxRichness, REML = FALSE)

# Reduced model (remove e.g. speciesPool random effect)
mem_richness_noSP <- lmer(maxRichness ~ scenario * timeStep_z * maxZ_z +
                          (1|X) + (1|Y),
                          data = maxRichness, REML = FALSE)

anova(mem_richness, mem_richness_noSP)


# - Perform mixed-effects regression for maxShannon

maxShannon <- maxShannon %>%
  mutate(
    timeStep_z = as.numeric(scale(timeStep)),
    maxZ_z     = as.numeric(scale(maxZ))
  )

# MEM
mem_shannon_test <- lmer(maxDiversity ~ scenario * timeStep_z * maxZ_z +
                          (1|X) + (1|Y) + (1|speciesPool),
                          data = maxShannon)
summary(mem_shannon_test)
# R2
r2(mem_shannon_test)


