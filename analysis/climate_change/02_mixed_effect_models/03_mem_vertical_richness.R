# Fit mixed-effects model (MEM) for vertical species richness maximum value and position
# 1. Find minimum model structure with species pool and random forest as random effects
# 2. Model diagnostics (plots)

library(readr)
library(lme4)
library(nlme)
library(lmerTest)
library(dplyr)
library(merTools)
library(DHARMa)
library(lattice)
library(ggplot2)
library(patchwork)
library(glmmTMB)
library(lattice)
library(DHARMa)
library(MASS)
library(performance)
library(tidyr)
library(purrr)

base_dir <- file.path("../modve_data/modve_output/regua")
richness_path <- file.path(base_dir, "vertical_diversity_peak_cc_vs_no_cc.csv")
DirectoryPlots <- file.path("../modve_figs/climate_change/mem/vertical_richness")

dir.create(DirectoryPlots)

if(!file.exists(richness_path)) {
  stop(richness_path, " missing. Run analysis/climate_change/01_community/04_plot_vertical_diversity.py.")
}

# Load data
richness_peak <- read_csv(richness_path)

# Sort scenarios and arrange data
richness_peak <- richness_peak %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

richness_peak_filtered <- richness_peak %>%
  filter(Year >= 2020) %>%
  # Identify all Scenario-Year-ForestID-SpeciesPool combos where Richness < 2
  anti_join(
    richness_peak %>%
      filter(level_5 == "Richness", Count < 2) %>%
      dplyr::select(Scenario, Year, ForestID, SpeciesPool),
    by = c("Scenario", "Year", "ForestID", "SpeciesPool")
  )

richness_peak_filtered$Year_c <- (richness_peak_filtered$Year - mean(richness_peak_filtered$Year)) / sd(richness_peak_filtered$Year)

# Scenatio to
richness_peak_filtered <- richness_peak_filtered %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

# Check distribution of responses
pdf(file.path(DirectoryPlots, "richness_hist.pdf"))
hist(richness_peak_filtered[richness_peak_filtered["level_5"] == "Richness", "Count"]$Count)
dev.off() # -> Poisson

pdf(file.path(DirectoryPlots, "richness_pos_hist.pdf"))
hist(richness_peak_filtered[richness_peak_filtered["level_5"] == "Richness", "Height"]$Height)
dev.off() # -> Gaussian?

shapiro.test(richness_peak_filtered[richness_peak_filtered["level_5"] == "Richness", "Height"]$Height) # Non-normal -> use Gamma

# #################################################################################################

# - Richness
richness <- richness_peak_filtered %>% filter(level_5 == "Richness")

mem_richness <- glmmTMB(
  Count ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
  dispformula = ~ SpeciesPool,  # models variance structure
  data = richness,
  family = compois(link = "log"),
  REML = TRUE
)

# 1. check singular fit
performance::check_singularity(mem_richness)

# 2. Check convergence
mem_richness$sdr$pdHess  # Should be TRUE

summary(mem_richness)  # Look at random effect variances

# 3. Check gradient
mem_richness$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_richness$fit$convergence  # Should be 0

round(exp(fixef(mem_richness)$cond), 2)

# --- DIAGONSTICS

prefix <- "richness_peak_val_"

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_richness, n = 1000)

pdf(file.path(DirectoryPlots, paste0(prefix, "res_qq_dharma_compois.pdf")),
    width = 10, height = 5)
plot(simulationOutput)
dev.off() # --> no outlier issues when modelling dispersion

testDispersion(simulationOutput)

# -- CHeck residuals on different levels - heteroscedasticity across groups

mf <- model.frame(mem_richness)            # model frame used to fit mem_richness
df <- data.frame(
  resid   = resid(mem_richness, type = "pearson"),
  fitted  = fitted(mem_richness),
  Scenario = mf$Scenario,
  SpeciesPool = mf$SpeciesPool,
  ForestID = mf$ForestID
)

# Named vectors for clean relabeling
forest_labs <- c(
  "0" = "Forest 1",
  "1" = "Forest 2",
  "2" = "Forest 3"
)

scenario_labs <- c(
  "CC"    = "Climate change",
  "No CC" = "Baseline"
)

p1 <- ggplot(df, aes(x = resid, y = as.factor(SpeciesPool))) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(
    ~ ForestID,
    ncol = 3,
    labeller = labeller(ForestID = forest_labs)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Pearson residuals",
    y = "Species pool"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major.y = element_blank()
  )

p2 <- ggplot(df, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(
    ~ Scenario,
    ncol = 3,
    labeller = labeller(Scenario = scenario_labs)
  ) +
  labs(
    x = "Fitted values",
    y = "Pearson residuals"
  ) +
  theme_bw()

pdf(file.path(DirectoryPlots, "richness_res_grouped.pdf"), width = 7, height = 7)
print(p2 / p1)
dev.off()

# -- Normality of randm effects -> OK but very small sample size to make strong statement

# Extract random effects
species_pool <- ranef(mem_richness)$cond$SpeciesPool[,1]
species_pool_scenario <- ranef(mem_richness)$cond$SpeciesPool$ScenarioCC
forest_id <- ranef(mem_richness)$cond$ForestID[,1]

# Plot histograms
pdf(file.path(DirectoryPlots, paste0(prefix, "rf_sp_normal_glmmtmb.pdf")))
hist(species_pool)
dev.off()

pdf(file.path(DirectoryPlots, paste0(prefix, "rf_sp_cc_normal_glmmtmb.pdf")))
hist(species_pool_scenario)
dev.off()

pdf(file.path(DirectoryPlots, paste0(prefix, "diag_rf_forest_normal_glmmtmb.pdf")))
hist(forest_id)
dev.off()

# Test for normality
shapiro.test(species_pool) # -> Normal
shapiro.test(species_pool_scenario) # -> Normal
shapiro.test(forest_id) # -> too few values

# ----- Variance partitioning ----- #

# Residual variance approx. from non-dispersion model
meme_no_disp <- update(mem_richness, dispformula = ~1)

performance::r2_nakagawa(meme_no_disp)
icc <- performance::icc(meme_no_disp)
icc_unadj <- icc$ICC_unadjusted * 100

# Get variances
vc <- VarCorr(mem_richness)

# Get variances
var_species_intercept <- vc$cond$SpeciesPool["(Intercept)", "(Intercept)"]
var_species_slope     <- vc$cond$SpeciesPool["ScenarioCC", "ScenarioCC"]
var_forest             <- vc$cond$ForestID["(Intercept)", "(Intercept)"]

# Total random-effect variance
total_var <- var_species_intercept + var_species_slope + var_forest

# Compute proportions
prop_species_intercept <- (var_species_intercept / total_var) * icc_unadj
prop_species_slope     <- (var_species_slope / total_var) * icc_unadj
prop_forest            <- (var_forest / total_var) * icc_unadj

cat(
  paste0(
    "Sp: ", round(var_species_intercept, 2), "(", round(prop_species_intercept, 2), "%)\n",
    "Sp_slope: ", round(var_species_slope, 2), "(", round(prop_species_slope, 2), "%)\n",
    "Forest: ", round(var_forest, 2), "(", round(prop_forest, 2), "%)"
  )
)

# - Richness peak position

mem_richness_pos <- glmmTMB(
  Height ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
  dispformula = ~ Scenario + SpeciesPool,  # models variance structure
  data = richness,
  family = Gamma(link = "log"),
  REML = TRUE
)

# 0. check singular fit
performance::check_singularity(mem_richness_pos)

# Check convergence
mem_richness_pos$sdr$pdHess  # Should be TRUE

summary(mem_richness_pos)  # Look at random effect variances

# 3. Check gradient
mem_richness_pos$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_richness_pos$fit$convergence  # Should be 0

round(exp(fixef(mem_richness_pos)$cond), 2)

# Did the model improve at all? -> Yes better to estimate dispersions
mem_richness_pos_null <- update(mem_richness_pos, dispformula = ~1)
anova(mem_richness_pos_null, mem_richness_pos)

# --- DIAGONSTICS

prefix <- "richness_peak_pos_"

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_richness_pos, n = 1000)

pdf(file.path(DirectoryPlots, paste0(prefix, "res_qq_dharma_compois.pdf")),
    width = 10, height = 5)
plot(simulationOutput)
dev.off()

testDispersion(simulationOutput)

# -- CHeck residuals on different levels - heteroscedasticity across groups

mf <- model.frame(mem_richness_pos)            # model frame used to fit mem_richness
df <- data.frame(
  resid   = resid(mem_richness_pos, type = "pearson"),
  fitted  = fitted(mem_richness_pos),
  Scenario = mf$Scenario,
  SpeciesPool = mf$SpeciesPool,
  ForestID = mf$ForestID
)


# Named vectors for clean relabeling
forest_labs <- c(
  "0" = "Forest 1",
  "1" = "Forest 2",
  "2" = "Forest 3"
)

scenario_labs <- c(
  "CC"    = "Climate change",
  "No CC" = "Baseline"
)

p1 <- ggplot(df, aes(x = resid, y = as.factor(SpeciesPool))) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(
    ~ ForestID,
    ncol = 3,
    labeller = labeller(ForestID = forest_labs)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Pearson residuals",
    y = "Species pool"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major.y = element_blank()
  )

p2 <- ggplot(df, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(
    ~ Scenario,
    ncol = 3,
    labeller = labeller(Scenario = scenario_labs)
  ) +
  labs(
    x = "Fitted values",
    y = "Pearson residuals"
  ) +
  theme_bw()

pdf(file.path(DirectoryPlots, "richness_hgt_res_grouped.pdf"), width = 7, height = 7)
print(p2 / p1)
dev.off()

# -- Normality of randm effects -> OK but very small sample size to make strong statement

# Extract random effects
species_pool <- ranef(mem_richness_pos)$cond$SpeciesPool[,1]
species_pool_scenario <- ranef(mem_richness_pos)$cond$SpeciesPool$ScenarioCC
forest_id <- ranef(mem_richness_pos)$cond$ForestID[,1]

# Plot histograms
pdf(file.path(DirectoryPlots, paste0(prefix, "rf_sp_normal_glmmtmb.pdf")))
hist(species_pool)
dev.off()

pdf(file.path(DirectoryPlots, paste0(prefix, "rf_sp_cc_normal_glmmtmb.pdf")))
hist(species_pool_scenario)
dev.off()

pdf(file.path(DirectoryPlots, paste0(prefix, "diag_rf_forest_normal_glmmtmb.pdf")))
hist(forest_id)
dev.off()

# Test for normality
shapiro.test(species_pool) # -> Normal
shapiro.test(species_pool_scenario) # -> Normal
shapiro.test(forest_id) # -> too few values

# ----- Variance partitioning ----- #

# Residual variance approx. from non-dispersion model
meme_no_disp <- update(mem_richness_pos, dispformula = ~1)

r2 <- performance::r2_nakagawa(meme_no_disp)
icc <- performance::icc(meme_no_disp)
icc_unadj <- icc$ICC_unadjusted * 100

# Get variances
vc <- VarCorr(meme_no_disp)

# Get variances
var_species_intercept <- vc$cond$SpeciesPool["(Intercept)", "(Intercept)"]
var_species_slope     <- vc$cond$SpeciesPool["ScenarioCC", "ScenarioCC"]
var_forest             <- vc$cond$ForestID["(Intercept)", "(Intercept)"]

# Total random-effect variance
total_var <- var_species_intercept + var_species_slope + var_forest

# Compute proportions
prop_species_intercept <- (var_species_intercept / total_var) * icc_unadj
prop_species_slope     <- (var_species_slope / total_var) * icc_unadj
prop_forest            <- (var_forest / total_var) * icc_unadj

cat(
  paste0(
    "Conditional R2: ", round(r2$R2_conditional, 3), "\n",
    "Marginal R2: ", round(r2$R2_marginal, 3), "\n",
    "ICC unadjusted: ", round(icc_unadj, 3), "%\n",
    "Sp: ", round(var_species_intercept, 2), "(", round(prop_species_intercept, 2), "%)\n",
    "Sp_slope: ", round(var_species_slope, 2), "(", round(prop_species_slope, 2), "%)\n",
    "Forest: ", round(var_forest, 2), "(", round(prop_forest, 2), "%)"
  )
)