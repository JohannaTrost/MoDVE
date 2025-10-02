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

# mem_pos <- glmer(
#   Position ~ Scenario * Year_c +
#     (1 | ForestID) +
#     (Scenario + Year_c | SpeciesPool),
#   data = species_distr_stats,
#   family = Gamma(link = "log")
# )
#
# summary(mem_pos)
#
# # Diagnostic plot
# pdf(file.path(DirectoryPlots, "Diagnostics", "diag_pos_res_qq_dharma_v6.pdf"))
# library(DHARMa)
#
# # Create scaled residuals
# simulationOutput <- simulateResiduals(fittedModel = mem_pos, n = 1000)
# plot(simulationOutput)
# dev.off()

# -> Issue with model residuals



# --- Alternative

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/Diagnostics")

library(glmmTMB)

mem_pos <- glmmTMB(
  Position ~ Scenario * Year_c + (Scenario + Year_c | SpeciesPool) + (1 | ForestID),
  dispformula = ~ SpeciesPool + Scenario,  # models variance structure
  data = species_distr_stats,
  family = Gamma(link = "log"),
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
pdf(file.path(DirectoryPlots, "diag_pos_res_scenario_frand_v6.pdf"))
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

pdf(file.path(DirectoryPlots, "diag_pos_res_sp_frand_v6_byScenario.pdf"))

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

pdf(file.path(DirectoryPlots, "diag_pos_res_qq_dharma_glmmtmb_v6.pdf"))
library(DHARMa)

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_pos, n = 1000)

# Or full diagnostic plot (includes QQ plot + more)
plot(simulationOutput)
dev.off()