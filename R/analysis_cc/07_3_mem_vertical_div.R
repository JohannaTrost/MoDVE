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

# Scenatio to
div_peak_filtered <- div_peak_filtered %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

# Check distribution of responses
pdf(file.path(DirectoryPlots, "richness_hist.pdf"))
hist(div_peak_filtered[div_peak_filtered["level_5"] == "Richness", "Count"]$Count)
dev.off() # -> Poisson

pdf(file.path(DirectoryPlots, "richness_pos_hist.pdf"))
hist(div_peak_filtered[div_peak_filtered["level_5"] == "Richness", "Height"]$Height)
dev.off() # -> Gaussian?

shapiro.test(div_peak_filtered[div_peak_filtered["level_5"] == "Richness", "Height"]$Height) # Non-normal -> use Gamma

pdf(file.path(DirectoryPlots, "abundance_hist.pdf"))
hist(div_peak_filtered[div_peak_filtered["level_5"] == "Abundance", "Count"]$Count)
dev.off() # -> Poisson

pdf(file.path(DirectoryPlots, "abundance_pos_hist.pdf"))
hist(div_peak_filtered[div_peak_filtered["level_5"] == "Abundance", "Height"]$Height)
dev.off() # -> gamma
shapiro.test(div_peak_filtered[div_peak_filtered["level_5"] == "Abundance", "Height"]$Height)

pdf(file.path(DirectoryPlots, "shannon_hist.pdf"))
hist(div_peak_filtered[div_peak_filtered["level_5"] == "Shannon", "Count"]$Count)
dev.off()

shapiro.test(div_peak_filtered[div_peak_filtered["level_5"] == "Shannon", "Count"]$Count)
# -> Gamma

# #################################################################################################

# - Richness

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/VerticalDiversity/Diagnostics")
dir.create(DirectoryPlots)

richness <- div_peak_filtered %>% filter(level_5 == "Richness")

mem_richness <- glmmTMB(
  Count ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
  dispformula = ~ SpeciesPool,  # models variance structure
  data = richness,
  family = compois(link = "log"),
  REML = TRUE
)

# 0. check singular fit
performance::check_singularity(mem_richness)

# Check convergence
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

# - Shannon

alpha_div <- div_peak_filtered %>% filter(level_5 == "Shannon")

mem_alpha_div <- glmmTMB(
  Count ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
   dispformula = ~ SpeciesPool,  # models variance structure
  data = alpha_div,
  #family = Gamma(link = "log"),
  family = tweedie(link = "log"),
  REML = TRUE
)

# 0. check singular fit
performance::check_singularity(mem_alpha_div)

# Check convergence
mem_alpha_div$sdr$pdHess  # Should be TRUE

summary(mem_alpha_div)  # Look at random effect variances

# 3. Check gradient
mem_alpha_div$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_alpha_div$fit$convergence  # Should be 0

# Get confidence intervals and estimates on the link scale
ci_link <- confint(mem_alpha_div, parm = "beta_", method = "wald")
ci_link_cond <- ci_link[grepl("cond", rownames(ci_link)), ]
ci_response <- exp(ci_link_cond)

estimates_with_ci <- data.frame(
  Parameter = names(fixef(mem_alpha_div)$cond),
  Estimate = round(exp(fixef(mem_alpha_div)$cond), 2),
  Lower_CI = round(ci_response[, 1], 2),
  Upper_CI = round(ci_response[, 2], 2)
)
print(estimates_with_ci)

# --- DIAGONSTICS

prefix <- "alpha_div_peak_val_"

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_alpha_div, n = 1000)

pdf(file.path(DirectoryPlots, paste0(prefix, "res_qq_dharma_tweedie_id.pdf")),
    width = 10, height = 5)
plot(simulationOutput)
dev.off() # --> no outlier issues when modelling dispersion

testDispersion(simulationOutput)
testUniformity(simulationOutput)

# -- CHeck residuals on different levels - heteroscedasticity across groups

mf <- model.frame(mem_alpha_div)            # model frame used to fit mem_richness
df <- data.frame(
  resid   = resid(mem_alpha_div, type = "pearson"),
  fitted  = fitted(mem_alpha_div),
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

pdf(file.path(DirectoryPlots, "alpha_div_res_grouped.pdf"), width = 7, height = 7)
print(p2 / p1)
dev.off()

# -- Normality of randm effects -> OK but very small sample size to make strong statement

# Extract random effects
species_pool <- ranef(mem_alpha_div)$cond$SpeciesPool[,1]
species_pool_scenario <- ranef(mem_alpha_div)$cond$SpeciesPool$ScenarioCC
forest_id <- ranef(mem_alpha_div)$cond$ForestID[,1]

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
meme_no_disp <- update(mem_alpha_div, dispformula = ~1)

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
    "ICC unadjusted: ", round((icc$ICC_unadjusted / r2$R2_conditional) * 100, 2), "%\n",
    "Sp: ", round(var_species_intercept, 2), "(", round(prop_species_intercept, 2), "%)\n",
    "Sp_slope: ", round(var_species_slope, 2), "(", round(prop_species_slope, 2), "%)\n",
    "Forest: ", round(var_forest, 2), "(", round(prop_forest, 2), "%)"
  )
)

# - Shannon peak position

mem_alpha_div_pos <- glmmTMB(
  Height ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
  dispformula = ~ SpeciesPool,  # models variance structure
  data = alpha_div,
  family = tweedie(link = "log"),
  REML = TRUE
)

# 0. check singular fit
performance::check_singularity(mem_alpha_div_pos)

# Check convergence
mem_alpha_div_pos$sdr$pdHess  # Should be TRUE

summary(mem_alpha_div_pos)  # Look at random effect variances

# 3. Check gradient
mem_alpha_div_pos$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_alpha_div_pos$fit$convergence  # Should be 0

# Get confidence intervals and estimates on the link scale
ci_link <- confint(mem_alpha_div_pos, parm = "beta_", method = "wald")
ci_link_cond <- ci_link[grepl("cond", rownames(ci_link)), ]
ci_response <- exp(ci_link_cond)

estimates_with_ci <- data.frame(
  Parameter = names(fixef(mem_alpha_div_pos)$cond),
  Estimate = round(exp(fixef(mem_alpha_div_pos)$cond), 2),
  Lower_CI = round(ci_response[, 1], 2),
  Upper_CI = round(ci_response[, 2], 2)
)
print(estimates_with_ci)

# Did the model improve at all? -> Yes better to estimate dispersions
mem_alpha_div_pos_null <- update(mem_alpha_div_pos, dispformula = ~1)
mem_cc <- update(mem_alpha_div_pos, dispformula = ~ Scenario)
mem_sp_cc <- update(mem_alpha_div_pos, dispformula = ~ SpeciesPool + Scenario)
mem_sp <- update(mem_alpha_div_pos, dispformula = ~ SpeciesPool)
anova(mem_alpha_div_pos_null, mem_cc, mem_sp_cc, mem_sp)

# --- DIAGONSTICS

prefix <- "alpha_div_peak_pos_"

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_alpha_div_pos, n = 1000)

pdf(file.path(DirectoryPlots, paste0(prefix, "res_qq_dharma_tweedie_v2.pdf")),
    width = 10, height = 5)
plot(simulationOutput)
dev.off() # --> no outlier issues when modelling dispersion

testDispersion(simulationOutput)
testOutliers(simulationOutput)

# -- CHeck residuals on different levels - heteroscedasticity across groups

mf <- model.frame(mem_alpha_div_pos)            # model frame used to fit mem_richness
df <- data.frame(
  resid   = resid(mem_alpha_div_pos, type = "pearson"),
  fitted  = fitted(mem_alpha_div_pos),
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

pdf(file.path(DirectoryPlots, "alpha_div_hgt_res_grouped.pdf"), width = 7, height = 7)
print(p2 / p1)
dev.off()

# -- Normality of randm effects -> OK but very small sample size to make strong statement

# Extract random effects
species_pool <- ranef(mem_alpha_div_pos)$cond$SpeciesPool[,1]
species_pool_scenario <- ranef(mem_alpha_div_pos)$cond$SpeciesPool$ScenarioCC
forest_id <- ranef(mem_alpha_div_pos)$cond$ForestID[,1]

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
meme_no_disp <- update(mem_alpha_div_pos, dispformula = ~1)

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
    "ICC unadjusted: ", (icc$ICC_unadjusted / r2$R2_conditional) * 100, "%\n",
    "Sp: ", round(var_species_intercept, 2), "(", round(prop_species_intercept, 2), "%)\n",
    "Sp_slope: ", round(var_species_slope, 2), "(", round(prop_species_slope, 2), "%)\n",
    "Forest: ", round(var_forest, 2), "(", round(prop_forest, 2), "%)"
  )
)

# - Abundance peak

abundance <- div_peak_filtered %>% filter(level_5 == "Abundance")

abundance$Year_log <- log(abundance$Year_c)
abundance$Year_c2 <- abundance$Year_c^2

mem_abundance <- glmmTMB(
  Count ~ Scenario * Year_c2 + (Scenario + Year_c2 | SpeciesPool) + (Scenario | ForestID),
  #dispformula = ~ SpeciesPool + Scenario,  # models variance structure
  data = abundance,
  family = poisson(link = "log"),
  REML = TRUE
)

# 0. check singular fit
performance::check_singularity(mem_abundance)

# Check convergence
mem_abundance$sdr$pdHess  # Should be TRUE

summary(mem_abundance)  # Look at random effect variances

# 3. Check gradient
mem_abundance$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_abundance$fit$convergence  # Should be 0

# Get confidence intervals and estimates on the link scale
ci_link <- confint(mem_abundance, parm = "beta_", method = "wald")
ci_response <- exp(ci_link[,c(1, 2)])

estimates_with_ci <- data.frame(
  Parameter = names(fixef(mem_abundance)$cond),
  Estimate = round(exp(fixef(mem_abundance)$cond), 2),
  Lower_CI = round(ci_response[, 1], 2),
  Upper_CI = round(ci_response[, 2], 2)
)
print(estimates_with_ci)

# --- DIAGONSTICS

prefix <- "abundance_peak_val_"

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_abundance, n = 1000)

pdf(file.path(DirectoryPlots, paste0(prefix, "res_qq_dharma_pois_v3.pdf")),
    width = 10, height = 5)
plot(simulationOutput)
dev.off()

testDispersion(simulationOutput)
testDispersion(simulationOutput)

# -- CHeck residuals on different levels - heteroscedasticity across groups

mf <- model.frame(mem_abundance)            # model frame used to fit mem_richness
df <- data.frame(
  resid   = resid(mem_abundance, type = "pearson"),
  fitted  = fitted(mem_abundance),
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

pdf(file.path(DirectoryPlots, "abundance_res_grouped.pdf"), width = 7, height = 7)
print(p2 / p1)
dev.off()

# Did the model improve at all? -> Does not help! stick with non dispersion model
mem_abundance_null <- update(mem_abundance, dispformula = ~1)
mem_cc <- update(mem_abundance, dispformula = ~ Scenario)
mem_sp_cc <- update(mem_abundance, dispformula = ~ SpeciesPool + Scenario)
anova(mem_abundance_null, mem_cc, mem_sp_cc, mem_abundance)

# -- Normality of randm effects -> OK but very small sample size to make strong statement

# Extract random effects
species_pool <- ranef(mem_abundance)$cond$SpeciesPool[,1]
species_pool_scenario <- ranef(mem_abundance)$cond$SpeciesPool$ScenarioCC
forest_id <- ranef(mem_abundance)$cond$ForestID[,1]

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

r2 <- performance::r2_nakagawa(mem_abundance)
icc <- performance::icc(mem_abundance)
icc_unadj <- icc$ICC_unadjusted * 100

# Get variances
vc <- VarCorr(mem_abundance)

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
    "ICC unadjusted: ", round((icc$ICC_unadjusted / r2$R2_conditional) * 100, 2), "%\n",
    "Sp: ", round(var_species_intercept, 2), "(", round(prop_species_intercept, 2), "%)\n",
    "Sp_slope: ", round(var_species_slope, 2), "(", round(prop_species_slope, 2), "%)\n",
    "Forest: ", round(var_forest, 2), "(", round(prop_forest, 2), "%)"
  )
)

# - Abundance peak position

abundance_no_outliers <- abundance[-outls, ]

mem_abundance_pos <- glmmTMB(
  Height ~ Scenario * Year_c + (Scenario | SpeciesPool) + (Scenario | ForestID),
  dispformula = ~ SpeciesPool + Scenario,  # models variance structure
  data = abundance,
  family = tweedie(link = "log"),
  REML = TRUE
)

mem_abundance_pos_null <- update(mem_abundance_pos, dispformula = ~1)
mem_cc <- update(mem_abundance_pos, dispformula = ~ Scenario)
mem_sp <- update(mem_abundance_pos, dispformula = ~ SpeciesPool)
mem_sp_cc <- update(mem_abundance_pos, dispformula = ~ SpeciesPool + Scenario)
anova(mem_abundance_pos_null, mem_cc, mem_sp, mem_sp_cc)

# 0. check singular fit
performance::check_singularity(mem_abundance_pos)

# Check convergence
mem_abundance_pos$sdr$pdHess  # Should be TRUE

summary(mem_abundance_pos)  # Look at random effect variances

# 3. Check gradient
mem_abundance_pos$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_abundance_pos$fit$convergence  # Should be 0

round(exp(fixef(mem_abundance_pos)$cond), 2)

# Get 95% confidence intervals on the link (log) scale
ci_link <- confint(mem_abundance_pos, parm = "beta_", level = 0.95)
ci_resp <- exp(ci_link)
round(as.data.frame(ci_resp), 2)

# --- DIAGONSTICS

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/VerticalDiversity/Diagnostics")
prefix <- "abundance_peak_pos_"

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_abundance_pos, n = 1000)

pdf(file.path(DirectoryPlots, paste0(prefix, "res_qq_dharma_tweedie.pdf")),
    width = 10, height = 5)
plot(simulationOutput)
dev.off()

testDispersion(simulationOutput)

# Quantify outliers
# outls_2 <- outliers(simulationOutput)
#
# abundance[outls_2, ] %>%
#   count(ForestID, sort = TRUE)
#
# abundance[outls_2, ] %>%
#   count(SpeciesPool, sort = TRUE)
#
# abundance[outls_2, ] %>%
#   count(Scenario, sort = TRUE)

# -- CHeck residuals on different levels - heteroscedasticity across groups

mf <- model.frame(mem_abundance_pos)            # model frame used to fit mem_richness
df <- data.frame(
  resid   = resid(mem_abundance_pos, type = "pearson"),
  fitted  = fitted(mem_abundance_pos),
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

pdf(file.path(DirectoryPlots, "abundance_hgt_res_grouped.pdf"), width = 7, height = 7)
print(p2 / p1)
dev.off()

# ----- Variance partitioning ----- #

# Residual variance approx. from non-dispersion model
meme_no_disp <- update(mem_abundance_pos, dispformula = ~1)

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
    "ICC unadjusted: ", round((icc$ICC_unadjusted / r2$R2_conditional) * 100, 2), "%\n",
    "Sp: ", round(var_species_intercept, 2), "(", round(prop_species_intercept, 2), "%)\n",
    "Sp_slope: ", round(var_species_slope, 2), "(", round(prop_species_slope, 2), "%)\n",
    "Forest: ", round(var_forest, 2), "(", round(prop_forest, 2), "%)"
  )
)