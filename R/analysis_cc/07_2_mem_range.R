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

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/Range")
dir.create(DirectoryPlots)

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

# Remove all where range is 0
species_distr_stats <- species_distr_stats %>%
  filter(Range > 0 & Year >= 2020)

pdf(file.path(DirectoryPlots, "diag_range_hist.pdf"))
hist(species_distr_stats$Range)
dev.off()

pdf(file.path(DirectoryPlots, "diag_iqr_hist.pdf"))
hist(species_distr_stats$IQR)
dev.off()

# #################################################################################################
#                               Mixed effects model for range                                  #
# #################################################################################################

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/Range/Diagnostics")
dir.create(DirectoryPlots)

mem_range <- glmer(
  Range ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
  #dispformula = ~ SpeciesPool + Scenario,  # models variance structure
  data = species_distr_stats,
  family = Gamma(link = "log"),
)

# check singular fit
isSingular(mem_range)

summary(mem_range)  # Look at random effect variances

# --- DIAGONSTICS

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_range, n = 1000)

pdf(file.path(DirectoryPlots, "diag_range_res_qq_dharma_v1.pdf"),
    width = 10, height = 5)
plot(simulationOutput)
dev.off()

# - Check heteroscedasrticity

# Get unique Scenarios
scenarios <- unique(species_distr_stats$Scenario)
levs <- levels(as.factor(species_distr_stats$SpeciesPool))
x <- as.numeric(species_distr_stats$SpeciesPool)
res <- residuals(mem_range, type = "pearson")

pdf(file.path(DirectoryPlots, "diag_range_res_byScenario_v1.pdf"),
    width = 5, height = 5)
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

# -- CHeck residuals on different levels

mf <- model.frame(mem_range)            # model frame used to fit mem_range
df <- data.frame(
  resid   = resid(mem_range, type = "pearson"),
  fitted  = fitted(mem_range),
  Scenario = mf$Scenario,
  SpeciesPool = mf$SpeciesPool,
  ForestID = mf$ForestID
)

# Scenario level plots
pdf(file.path(DirectoryPlots, "diag_range_res_spread_scenario_v1.pdf"))
xyplot(resid ~ fitted | Scenario, data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# Species pool level plots
pdf(file.path(DirectoryPlots, "diag_range_res_spread_sp_v1.pdf"))
xyplot(resid ~ fitted | as.factor(SpeciesPool), data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# Forest level plots
pdf(file.path(DirectoryPlots, "diag_range_res_spread_forest_v1.pdf"))
xyplot(resid ~ fitted | as.factor(ForestID), data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# ---- Account for heteroscedasticity in SP

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/Range/Diagnostics/v2_with_dispformula")
dir.create(DirectoryPlots)

mem_range <- glmmTMB(
  Range ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
  dispformula = ~ Scenario + SpeciesPool,  # models variance structure
  data = species_distr_stats,
  family = Gamma(link = "log"),
  REML = TRUE
)

# 0. check singular fit
performance::check_singularity(mem_range)

# Check convergence
mem_range$sdr$pdHess  # Should be TRUE

summary(mem_range)  # Look at random effect variances

# 3. Check gradient
mem_range$sdr$gradient.fixed  # Should be close to zero

# 4. Convergence code
mem_range$fit$convergence  # Should be 0

# --- DIAGONSTICS

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_range, n = 1000)

pdf(file.path(DirectoryPlots, "diag_range_res_qq_dharma.pdf"),
    width = 10, height = 5)
plot(simulationOutput)
dev.off()


# -- CHeck residuals on different levels

mf <- model.frame(mem_range)            # model frame used to fit mem_range
df <- data.frame(
  resid   = resid(mem_range, type = "pearson"),
  fitted  = fitted(mem_range),
  Scenario = mf$Scenario,
  SpeciesPool = mf$SpeciesPool,
  ForestID = mf$ForestID
)

# Scenario level plots
pdf(file.path(DirectoryPlots, "diag_range_res_spread_scenario_v1.pdf"))
xyplot(resid ~ fitted | Scenario, data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# Species pool level plots
pdf(file.path(DirectoryPlots, "diag_range_res_spread_sp_v1.pdf"))
xyplot(resid ~ fitted | as.factor(SpeciesPool), data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# Forest level plots
pdf(file.path(DirectoryPlots, "diag_range_res_spread_forest_v1.pdf"))
xyplot(resid ~ fitted | as.factor(ForestID), data = df,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(h = 0)   # horizontal line at 0
       },
       xlab = "Fitted values", ylab = "Pearson residuals")
dev.off()

# Did the model improve at all?

mem_range_null <- update(mem_range, dispformula = ~1)
anova(mem_range_null, mem_range)
