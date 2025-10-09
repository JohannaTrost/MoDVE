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
species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc")

# Get stats
species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
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

pdf(file.path(DirectoryPlots, "diag_pos_hist.pdf"))
hist(species_distr_stats$Position)
dev.off()

# #################################################################################################
#                               Mixed effects model for position                                  #
# #################################################################################################

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/Diagnostics")

mem_pos <- glmmTMB(
  Position ~ Scenario * Year_c + (1 | SpeciesPool) + (1 | ForestID),
  #dispformula = ~ SpeciesPool + Scenario,  # models variance structure
  data = species_distr_stats,
  family = Gamma(link = "log"),
  REML = TRUE
)

# 0. check singular fit
performance::check_singularity(mem_pos)

# 1. Check convergence
mem_pos$sdr$pdHess  # Should be TRUE

# 2. Check for singular fit (near-zero variances)
summary(mem_pos)  # Look at random effect variances
VarCorr(mem_pos)  # Detailed variance-covariance

# 3. Check gradient
mem_pos$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_pos$fit$convergence  # Should be 0

# Extract scaled residuals
res <- residuals(mem_pos, type = "pearson")
levs <- levels(species_distr_stats$Scenario)
x <- as.numeric(species_distr_stats$Scenario)

# Make the plot
pdf(file.path(DirectoryPlots, "diag_pos_res_scenario_frand_glmmtmb_v6.pdf"))
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

pdf(file.path(DirectoryPlots, "diag_pos_res_sp_frand_glmmtmb_v6_byScenario.pdf"))

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

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_pos, n = 1000)

pdf(file.path(DirectoryPlots, "diag_pos_res_qq_dharma_glmmtmb_v6.pdf"))
plot(simulationOutput)
dev.off()

# -> issues

# -- CHeck residuals on different levels

mf <- model.frame(mem_pos)            # model frame used to fit mem_pos
df <- data.frame(
  resid   = resid(mem_pos, type = "pearson"),
  fitted  = fitted(mem_pos),
  Scenario = mf$Scenario,
  SpeciesPool = mf$SpeciesPool,
  ForestID = mf$ForestID
)

# Scenario level plots
pdf(file.path(DirectoryPlots, "diag_pos_res_spread_scenario_glmmtmb_v6.pdf"))

xyplot(resid ~ fitted | Scenario, data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# Species pool level plots
pdf(file.path(DirectoryPlots, "diag_pos_res_spread_sp_glmmtmb_v6.pdf"))
xyplot(resid ~ fitted | as.factor(SpeciesPool), data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# Forest level plots
pdf(file.path(DirectoryPlots, "diag_pos_res_spread_forest_glmmtmb_v6.pdf"))
xyplot(resid ~ fitted | as.factor(ForestID), data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# -- Normality of randm effects

# Extract random effects
species_pool <- ranef(mem_pos)$cond$SpeciesPool[,1]
forest_id <- ranef(mem_pos)$cond$ForestID[,1]

# Plot histograms
pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_normal_glmmtmb.pdf"))
hist(species_pool)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_forest_normal_glmmtmb.pdf"))
hist(forest_id)
dev.off()

# Test for normality
shapiro.test(species_pool)
shapiro.test(forest_id)

# ----- Variance partitioning -----

VarCorr(mem_pos)

performance::r2_nakagawa(mem_pos)

# Try glmer

mem_pos_v2 <- glmer(Position ~ Scenario * Year_c + (1 | SpeciesPool) + (1 | ForestID),
                    family = Gamma(link = log),
                    data = species_distr_stats)

# Check sinular fit
performance::check_singularity(mem_pos_v2)

# Check multicpllineratity issue -> values reasonably low -> no predictor issues
car::vif(lm(Position ~ Scenario * Year_c, data=species_distr_stats))

summary(mem_pos_v2)

# Variance partitioning

icc(mem_pos_v2) # -> random effects explain <1% of variance

# Variance components
sp_var <- 0.0003057
frst_var <- 0.0009625
total_var <- sp_var + frst_var + sigma(mem_pos_v2)^2
sp_vc <- 100*(sp_var / total_var)
frst_vc <- 100*(frst_var / total_var)

print(paste0("Species pool variance component: ", round(sp_vc, 4)))
print(paste0("Forest ID variance component: ", round(frst_vc, 4)))
