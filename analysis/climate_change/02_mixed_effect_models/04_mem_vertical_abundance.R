# Fit mixed-effects model (MEM) for vertical species abundance maximum value and position
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
abundance_path <- file.path(base_dir, "vertical_diversity_peak_cc_vs_no_cc.csv")
DirectoryPlots <- file.path("../modve_figs/climate_change/mem/vertical_abundance")

dir.create(DirectoryPlots)

if(!file.exists(abundance_path)) {
  stop(abundance_path, " missing. Run analysis/climate_change/01_community/04_plot_vertical_diversity.py.")
}

# Load data
abundance_peak <- read_csv(abundance_path)

# Sort scenarios and arrange data
abundance_peak <- abundance_peak %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

abundance_peak_filtered <- abundance_peak %>%
  filter(Year >= 2020) %>%
  # Identify all Scenario-Year-ForestID-SpeciesPool combos where Richness < 2
  anti_join(
    abundance_peak %>%
      filter(level_5 == "Richness", Count < 2) %>%
      dplyr::select(Scenario, Year, ForestID, SpeciesPool),
    by = c("Scenario", "Year", "ForestID", "SpeciesPool")
  )

abundance_peak_filtered$Year_c <- (abundance_peak_filtered$Year - mean(abundance_peak_filtered$Year)) / sd(abundance_peak_filtered$Year)

# Scenatio to
abundance_peak_filtered <- abundance_peak_filtered %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

pdf(file.path(DirectoryPlots, "abundance_hist.pdf"))
hist(abundance_peak_filtered[abundance_peak_filtered["level_5"] == "Abundance", "Count"]$Count)
dev.off() # -> Poisson

pdf(file.path(DirectoryPlots, "abundance_pos_hist.pdf"))
hist(abundance_peak_filtered[abundance_peak_filtered["level_5"] == "Abundance", "Height"]$Height)
dev.off() # -> gamma
shapiro.test(abundance_peak_filtered[abundance_peak_filtered["level_5"] == "Abundance", "Height"]$Height)

# - Abundance peak

abundance <- abundance_peak_filtered %>% filter(level_5 == "Abundance")

abundance$Year_log <- log(abundance$Year_c)
abundance$Year_c2 <- abundance$Year_c^2

mem_abundance <- glmmTMB(
  Count ~ Scenario * Year_c2 + (Scenario + Year_c2 | SpeciesPool) + (Scenario | ForestID),
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