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

# Repsonse vars:
# - Avg. species position
# - Species range size
# - Species richness (peak)
# - Species abundance (peak)
# - Species diversity (peak)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")
species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

# Get stats
species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, ForestID, TimeStep, Year, ) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop")

# Scale Year
species_distr_stats$Year_c <- (species_distr_stats$Year - mean(species_distr_stats$Year)) / sd(species_distr_stats$Year)

# Sort scenarios and arrange data
species_distr_stats <- species_distr_stats %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

mem_pos <- lmer(Position ~ Scenario * Year_c + (1 | ForestID) + (Scenario | SpeciesPool),
                data = species_distr_stats)

# Removed random slope for forestID due to convergence issues and
# almost perfect correlation of random slope and intercept

summary(mem_pos)

# > summary(mem_pos)
# Linear mixed model fit by REML ['lmerMod']
# Formula: Position ~ Scenario * Year_c + (1 | ForestID) + (Scenario | SpeciesPool)
#    Data: species_distr_stats
#
# REML criterion at convergence: 33292.9
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -3.3445 -0.5968 -0.0885  0.5604  6.0293
#
# Random effects:
#  Groups      Name        Variance Std.Dev. Corr
#  SpeciesPool (Intercept) 2.036    1.427            -> moderate variation of position among species pools
#              ScenarioCC  1.162    1.078    0.11    -> - corr: species pools starting higher are only slightly more
#                                                               likely to experience stronger downward shifts under climate change
#                                                       - Species pools sensitivity to CC varies, effect of CC varies by +-1m
#                                                       (interpretation: different combinations of tolerances translate into differential
#                                                                        vertical responses to climate change)
#  ForestID    (Intercept) 1.295    1.138            -> some variation of position among forests
#  Residual                5.922    2.434
# Number of obs: 7185, groups:  SpeciesPool, 10; ForestID, 3
#
# Fixed effects:
#                    Estimate Std. Error t value
# (Intercept)       19.905605   0.798142  24.940 -> Occur height at 19.9 m year=0 and no CC
# ScenarioCC        -0.610442   0.345734  -1.766 -> Occur lower (0.61 m), but marginal effect (t-value < 2)
# Year_c            -0.053581   0.001171 -45.761 -> Strong decrease over time per unit of centered year
# ScenarioCC:Year_c -0.024469   0.001661 -14.733 -> With climate change they additional move 0.024m lower per year
#
# Correlation of Fixed Effects:
#             (Intr) ScnrCC Year_c
# ScenarioCC   0.054
# Year_c       0.000  0.000
# ScnrCC:Yr_c  0.000  0.000 -0.705

round(vcov(mem_pos), 2)
VarCorr(mem_pos)

# -- Is the CC effect significant?

# 1
mem_pos <- lmer(Position ~ Scenario * Year_c + as.factor(ForestID) + (Scenario|SpeciesPool),
                data = species_distr_stats)
summary(mem_pos)
#                     Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)        1.991e+01  7.981e-01  4.130e+00  24.940 1.16e-05 ***
# ScenarioCC        -6.104e-01  3.457e-01  9.001e+00  -1.766    0.111
# Year_c            -5.358e-02  1.171e-03  7.161e+03 -45.761  < 2e-16 ***
# ScenarioCC:Year_c -2.447e-02  1.661e-03  7.161e+03 -14.733  < 2e-16 ***

# 2
confint(mem_pos, method = "Wald")  # or method = "profile" for more accurate intervals

# 3 nested models
model_full <- lmer(Position ~ Scenario + Year_c + (1|ForestID) + (1|SpeciesPool), data = species_distr_stats)
model_reduced <- lmer(Position ~ Year_c + (1|ForestID) + (1|SpeciesPool), data = species_distr_stats)
anova(model_reduced, model_full) # -> Significant effect of Scenario (< 2.2e-16 ***)

# ------------------ Predictions for Viz ------------------ #

nScenarios <- length(unique(species_distr_stats$Scenario))
nForests <- length(unique(species_distr_stats$ForestID))
nSpeciesPools <- length(unique(species_distr_stats$SpeciesPool))
nYears <- length(seq(1981, 2100))

newdata <- data.frame(
  Year = rep(seq(1981, 2100), nScenarios * nForests * nSpeciesPools),
  Scenario = c(rep("No CC", nYears * nForests * nSpeciesPools),
               rep("CC", nYears * nForests * nSpeciesPools)),
  ForestID = as.factor(rep(c(0, 1, 2), nYears * nSpeciesPools * nScenarios)),
  SpeciesPool = as.factor(rep(seq(10), nYears * nForests * nScenarios))
)
newdata$Year_c <- (newdata$Year - mean(species_distr_stats$Year)) / sd(species_distr_stats$Year)
newdata$Position <- predict(mem_pos, newdata)
# -> not working becaus new combinations of REs

# CI with uncertainty from random effects
# Prediction intervals (including RE uncertainty)
preds <- predictInterval(mem_pos,
                         newdata = species_distr_stats,
                         level = 0.95,
                         n.sims = 1000,     # number of draws
                         which = "full",    # includes RE uncertainty
                         include.resid.var = TRUE)

cbind(preds, species_distr_stats) %>%
  tibble(.) %>%
  mutate(Year_c = purrr::map_dbl(Year_c, 1)) %>%
  write_csv(., file.path(base_dir, "a5_pred_pos_cc_vs_no_cc.csv"))

# ------------------ DIAGNOSTICS ------------------ #

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/Diagnostics")

# 1. Residual checks

# res <- simulateResiduals(fittedModel = mem_pos, n = 1000)
#
# # Default diagnostic plots
# pdf(file.path(DirectoryPlots, "diag_pos_rf_dharma.pdf"))
# plot(res)
# dev.off()

# Additional checks
# testUniformity(res)   # Tests if residuals are uniformly distributed
# testDispersion(res)   # Tests for over/underdispersion
# testZeroInflation(res) # For count models

# 1. Normality of randm effects

# Extract random effects
species_pool <- ranef(mem_pos)$SpeciesPool[,1]
species_pool_slope <- ranef(mem_pos)$SpeciesPool[,2]
forest_id <- ranef(mem_pos)$ForestID[,1]

# Plot histograms
pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_normal.pdf"))
hist(species_pool)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_slope_normal.pdf"))
hist(species_pool_slope)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_forest_normal.pdf"))
hist(forest_id)
dev.off()

# Test for normality
shapiro.test(species_pool) # W = 0.9199, p-value = 0.3561 -> normal
shapiro.test(species_pool_slope) # W = 0.92966, p-value = 0.4446 -> normal
shapiro.test(forest_id) # W = 0.81578, p-value = 0.1527 -> normal

# 2. QQ plot of residuals
pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_qq.pdf"))
qqmath(species_pool,
       id = 0.05,
       panel = function(x, ...) {
         panel.qqmath(x, ...)                 # default Q-Q plot
         panel.abline(a = 0, b = 1, col = "black")  # 45° identity line
       })
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_slope_qq.pdf"))
qqmath(species_pool_slope,
       id = 0.05,
       panel = function(x, ...) {
         panel.qqmath(x, ...)                 # default Q-Q plot
         panel.abline(a = 0, b = 1, col = "black")  # 45° identity line
       })
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_forest_qq.pdf"))
qqmath(forest_id,
       id = 0.05,
       panel = function(x, ...) {
         panel.qqmath(x, ...)                 # default Q-Q plot
         panel.abline(a = 0, b = 1, col = "black")  # 45° identity line
       })
dev.off()

# 2. Check distribution of Residuals

pdf(file.path(DirectoryPlots, "diag_pos_res_overall.pdf"))
plot(mem_pos)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_res_qq.pdf"))
qqmath(mem_pos, id=0.05) # set significance level for outlier test
dev.off()

# - Plot residuals grouped by Scenario

# Extract scaled residuals
res <- resid(mem_pos, scaled = TRUE)
levs <- levels(species_distr_stats$Scenario)
x <- as.numeric(species_distr_stats$Scenario)

# Make the plot
pdf(file.path(DirectoryPlots, "diag_pos_res_scenario.pdf"))
plot(x, res,
     xaxt = "n",  # suppress default x-axis
     xlab = "Scenario",
     ylab = "Scaled Residuals",
     main = "Residuals vs Scenario")
abline(h = 0, col = "red", lty = 2)
axis(1, at = 1:length(levs), labels = levs)
dev.off()

# Residual spread by Scenario
pdf(file.path(DirectoryPlots, "diag_pos_res_spread_scenario.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | Scenario, abline = 0)
dev.off()
# -> this violates the assumption of homogeneous (constant) variance along a predictor.

pdf(file.path(DirectoryPlots, "diag_pos_res_spread_forest.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | ForestID, abline = 0)
dev.off()

# ------ Prepro ------ #

species_distr_stats$ForestID <- as.factor(species_distr_stats$ForestID)
species_distr_stats$SpeciesPool <- as.factor(species_distr_stats$SpeciesPool)
species_distr_stats %>% filter(Year > 2000) -> species_distr_stats

mem_pos_lme_var <- lme(Position ~ Scenario * Year_c,
                   random = list(SpeciesPool=~Scenario+Year_c,
                                 ForestID=~Scenario+Year_c),
                   #weights = varIdent(form = ~1 | Scenario),
                   data = species_distr_stats,
                   method = "REML")
# weights = varComb( varIdent(form = ~1 | Scenario),
#                                       varIdent(form = ~1 | Year_c))

# -> no convergence
# mem_pos_lme_max_var <- lme(Position ~ Scenario * Year_c * ForestID,
#                    random = ~Scenario + ForestID + Year_c | SpeciesPool,
#                    weights = varComb( varIdent(form = ~1 | Scenario),
#                                       varIdent(form = ~1 | ForestID)),
#                    data = species_distr_stats,
#                    method = "REML")
# -> Also with var structure for SpeciesPools it did not converfe

# Best model structure
mem_pos_lme <- lme(Position ~ Scenario * Year_c,
                   random = list(SpeciesPool=~Year_c,
                                 ForestID=~Scenario+Year_c),
                   #weights = varIdent(form = ~1 | Scenario),
                   data = species_distr_stats,
                   method = "REML")

# A more complex random structure did not converge
mem_pos <- lmer(Position ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
                data = species_distr_stats)

summary(mem_pos)
# - Main take away from summary
# Climate change effect: Species position is lower under climate change, and the gap widens over time.
# Context dependence: Forests differ strongly in how species respond, and these differences also play out over time.
# Random effects: Variation among species pools in both baseline position and climate sensitivity is meaningful, showing real heterogeneity.

# ------------------ 2. DIAGNOSTICS ------------------ #

# 1. - Heteroscedasticity across Scenarios and Forests

# Extract scaled residuals
res <- resid(mem_pos, scaled = TRUE)
levs <- levels(species_distr_stats$Scenario)
x <- as.numeric(species_distr_stats$Scenario)

# Make the plot
pdf(file.path(DirectoryPlots, "diag_pos_res_scenario_frand_v4.pdf"))
plot(x, res,
     xaxt = "n",  # suppress default x-axis
     xlab = "Scenario",
     ylab = "Scaled Residuals",
     main = "Residuals vs Scenario")
abline(h = 0, col = "red", lty = 2)
axis(1, at = 1:length(levs), labels = levs)
dev.off()

# - Residuals by SpeciesPool and Scenario
levs <- levels(as.factor(species_distr_stats$SpeciesPool))
x <- as.numeric(species_distr_stats$SpeciesPool)

pdf(file.path(DirectoryPlots, "diag_pos_res_sp_frand_v4_byScenario.pdf"))

# Get unique Scenarios
scenarios <- unique(species_distr_stats$Scenario)

par(mfrow = c(length(scenarios), 1), mar = c(6, 4, 3, 1))  # one row per scenario
for (sc in scenarios) {
  idx <- species_distr_stats$Scenario == sc
  plot(x[idx], res[idx],
       xaxt = "n",
       xlab = "Species pool",
       ylab = "Scaled Residuals",
       main = paste("Residuals vs species pools - Scenario:", sc))
  abline(h = 0, col = "red", lty = 2)
  axis(1, at = 1:length(levs), labels = levs, las = 2)
}

dev.off()

# Residual spread by Scenario
pdf(file.path(DirectoryPlots, "diag_pos_res_spread_scenario_frand_v4.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | Scenario, abline = 0)
dev.off()
# -> this violates the assumption of homogeneous (constant) variance along a predictor.

pdf(file.path(DirectoryPlots, "diag_pos_res_spread_forest_frand_v4.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | ForestID, abline = 0)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_res_spread_sp_frand_v4.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | as.factor(SpeciesPool), abline = 0)
dev.off()

# - VarCov RFs

#getVarCov(mem_pos, type = "random.effects")

# Get variance for each Scenario
#cc_var <- 0 # TODO
#cat("CC scenario variance estimate: ", mem_pos$sigma^2 * cc_var^2)

# - Relative heteroscedasticity: The CC scenario is noisier — your model suggests the spread of residuals is larger there than in the No CC scenario.
# - Implication: Predictions for the CC scenario will be less precise (larger error spread).
# - Biological/experimental interpretation: If CC is some experimental manipulation, it increases variability in responses, not just mean shifts.

# 2. - Normality of random effects

# Extract random effects
species_pool <- ranef(mem_pos)$SpeciesPool[,1]
sp_slope_scenario <- ranef(mem_pos)$SpeciesPool$ScenarioCC
sp_frst <- ranef(mem_pos)$ForestID[,1]

# Plot histograms
pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_normal_frand_v4.pdf"))
hist(species_pool)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_slope_scenario_normal_frand_v4.pdf"))
hist(sp_slope_scenario)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_frst_normal_frand_v4.pdf"))
hist(sp_frst)
dev.off()

# Test for normality
shapiro.test(species_pool) # W = 0.96094, p-value = 0.7966 -> normal
shapiro.test(sp_slope_scenario) # W = 0.92941, p-value = 0.4421 -> normal
shapiro.test(sp_frst) # W = 0.90383, p-value = 0.2413 -> normal

# 3. - QQ plot of random effects

# 1. QQ plot of residuals
pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_qq_frand_v4.pdf"))
qqmath(species_pool,
       id = 0.05,
       panel = function(x, ...) {
         panel.qqmath(x, ...)                 # default Q-Q plot
         panel.abline(a = 0, b = 1, col = "black")  # 45° identity line
       })
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_slope_scenario_qq_frand_v4.pdf"))
qqmath(sp_slope_scenario,
       id = 0.05,
       panel = function(x, ...) {
         panel.qqmath(x, ...)                 # default Q-Q plot
         panel.abline(a = 0, b = 1, col = "black")  # 45° identity line
       })
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_frst_qq_frand_v4.pdf"))
qqmath(sp_frst,
       id = 0.05,
       panel = function(x, ...) {
         panel.qqmath(x, ...)                 # default Q-Q plot
         panel.abline(a = 0, b = 1, col = "black")  # 45° identity line
       })
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_slope_frst2_qq_frand_v4.pdf"))
qqmath(sp_slope_frst2,
       id = 0.05,
       panel = function(x, ...) {
         panel.qqmath(x, ...)                 # default Q-Q plot
         panel.abline(a = 0, b = 1, col = "black")  # 45° identity line
       })
dev.off()

# 2. Check distribution of Residuals

pdf(file.path(DirectoryPlots, "diag_pos_res_overall_frand_v4.pdf"))
plot(mem_pos)
dev.off()

# QQ plot
pdf(file.path(DirectoryPlots, "diag_pos_res_qq_frand_v4.pdf"))

# Extract residuals
resids <- resid(mem_pos)

# Create Q-Q plot with line
qqnorm(resids, main = "Normal Q-Q Plot of Standardized Residuals")
qqline(resids, col = "red", lwd = 2)

# Identify outliers and annotate with SpeciesPool
outlier_threshold <- 3
outliers <- which(abs(scale(resids)) > outlier_threshold)

# Get theoretical quantiles for outliers
qq_data <- qqnorm(resids, plot.it = FALSE)

# Add text labels for SpeciesPool
text(qq_data$x[outliers],
     resids[outliers],
     labels = species_distr_stats$SpeciesPool[outliers],
     col = "darkgrey",
     cex = 0.7,
     pos = 4)  # pos = 4 places text to the right of points

dev.off()

# Investigate outliers

# Or identify programmatically
resid_std <- residuals(mem_pos, type = "pearson") /
             sd(residuals(mem_pos, type = "pearson"))
outliers <- which(abs(resid_std) > 3)

# Examine them
species_distr_stats[outliers, ]

# Are they from specific groups?
table(species_distr_stats$SpeciesPool[outliers])
table(species_distr_stats$Scenario[outliers])
# -> SpeciesPool 9 (and to a lesser extent 1) have unique temporal dynamics under climate change in later years that need to be modeled

# Consdier transforming Position
pdf(file.path(DirectoryPlots, "diag_pos_distr.pdf"))
hist(species_distr_stats$Position)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_distr.pdf"))
hist(species_distr_stats$log_Position)
dev.off()

shapiro.test(sample(species_distr_stats$log_Position, 5000))

# --- To overcome all hurdles (crossed + variance structure + 3 forests only) --- #


library(brms)

get_prior(Position ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID), data=species_distr_stats)

mem_pos <- brm(
  bf(Position ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
     sigma ~ Scenario),
  data = species_distr_stats, family = gaussian(),
  # Consider priors for ForestID to regularize the 3-level variance
  prior = c(
    prior(cauchy(0, 1), class = sd)
  ),
  chains = 4, iter = 5000, warmup = 1000,
  save_pars = save_pars(all=TRUE) # for marginal likelihood
)

# --- Alternative

library(glmmTMB)

mem_pos <- glmmTMB(
  Position ~ Scenario * Year_c + (Scenario + Year_c | SpeciesPool) + (1 | ForestID),
  dispformula = ~ SpeciesPool + Scenario,  # models variance structure
  data = species_distr_stats,
  REML = TRUE
)

# 1. Check convergence
mem_pos$sdr$pdHess  # Should be TRUE

# 2. Check for singular fit (near-zero variances)
summary(mem_pos)  # Look at random effect variances
VarCorr(mem_pos)  # Detailed variance-covariance

# 3. Check gradient
mem_pos$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_pos$fit$convergence  # Should be 0

resids <- residuals(mem_pos, type = "pearson")

# Extract scaled residuals
res <- residuals(mem_pos, type = "pearson")
levs <- levels(species_distr_stats$Scenario)
x <- as.numeric(species_distr_stats$Scenario)

# Make the plot
pdf(file.path(DirectoryPlots, "diag_pos_res_scenario_frand_v5.pdf"))
plot(x, res,
     xaxt = "n",  # suppress default x-axis
     xlab = "Scenario",
     ylab = "Scaled Residuals",
     main = "Residuals vs Scenario")
abline(h = 0, col = "red", lty = 2)
axis(1, at = 1:length(levs), labels = levs)
dev.off()

# - Residuals by SpeciesPool and Scenario
levs <- levels(as.factor(species_distr_stats$SpeciesPool))
x <- as.numeric(species_distr_stats$SpeciesPool)

pdf(file.path(DirectoryPlots, "diag_pos_res_sp_frand_v5_byScenario.pdf"))

# Get unique Scenarios
scenarios <- unique(species_distr_stats$Scenario)

par(mfrow = c(length(scenarios), 1), mar = c(6, 4, 3, 1))  # one row per scenario
for (sc in scenarios) {
  idx <- species_distr_stats$Scenario == sc
  plot(x[idx], res[idx],
       xaxt = "n",
       xlab = "Species pool",
       ylab = "Scaled Residuals",
       main = paste("Residuals vs species pools - Scenario:", sc))
  abline(h = 0, col = "red", lty = 2)
  axis(1, at = 1:length(levs), labels = levs, las = 2)
}

dev.off()

# Residual spread by Scenario
pdf(file.path(DirectoryPlots, "diag_pos_res_spread_scenario_frand_v5.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | Scenario, abline = 0)
dev.off()

# -> this violates the assumption of homogeneous (constant) variance along a predictor.

pdf(file.path(DirectoryPlots, "diag_pos_res_spread_forest_frand_v5.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | ForestID, abline = 0)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_res_spread_sp_frand_v5.pdf"))
plot(mem_pos, resid(.) ~ fitted(.) | as.factor(SpeciesPool), abline = 0)
dev.off()


# QQ plot
pdf(file.path(DirectoryPlots, "diag_pos_res_qq_dharma_glmmtmb2_v5.pdf"))

# Create Q-Q plot with line
qqnorm(resids, main = "Normal Q-Q Plot of Standardized Residuals")
qqline(resids, col = "red", lwd = 2)

# Identify outliers and annotate with SpeciesPool
outlier_threshold <- 3
outliers <- which(abs(scale(resids)) > outlier_threshold)

# Get theoretical quantiles for outliers
qq_data <- qqnorm(resids, plot.it = FALSE)

# Add text labels for SpeciesPool
text(qq_data$x[outliers],
     resids[outliers],
     labels = species_distr_stats$SpeciesPool[outliers],
     col = "darkgrey",
     cex = 0.7,
     pos = 4)  # pos = 4 places text to the right of points

dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_res_qq_dharma_glmmtmb_v5.pdf"))
library(DHARMa)

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_pos, n = 1000)

# QQ plot
plotQQunif(simulationOutput)

# Or full diagnostic plot (includes QQ plot + more)
plot(simulationOutput)
dev.off()