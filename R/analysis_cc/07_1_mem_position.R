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

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/Position/Diagnostics")

# Try glmer
mem_pos <- glmer(Position ~ Scenario * Year_c + (1 | SpeciesPool) + (1 | ForestID),
                    family = Gamma(link = log),
                    data = species_distr_stats)

# Check sinular fit
performance::check_singularity(mem_pos)

# Check multicpllineratity issue -> values reasonably low -> no predictor issues
car::vif(lm(Position ~ Scenario * Year_c, data=species_distr_stats))

summary(mem_pos)

# Estimates on response scal e
ci_link <- confint(mem_pos, parm = "beta_", method = "Wald")
# Exponentiate to get confidence intervals on the response scale
ci_response <- exp(ci_link)
# summary on response scale
results <- data.frame(
  Parameter = rownames(ci_response),
  Estimate = round(exp(fixef(mem_pos)), 2),
  Lower_CI = round(ci_response[, 1], 2),
  Upper_CI = round(ci_response[, 2], 2)
)
print(results, row.names = FALSE)

# Create scaled residuals
simulationOutput <- simulateResiduals(fittedModel = mem_pos, n = 1000)

pdf(file.path(DirectoryPlots, "diag_pos_res_qq_dharma_v2.pdf"),
    width = 10, height = 5)
plot(simulationOutput)
dev.off()

# -- CHeck residuals on different levels

mf <- model.frame(mem_pos)            # model frame used to fit mem_pos
df <- data.frame(
  resid   = resid(mem_pos, type = "pearson"),
  fitted  = fitted(mem_pos),
  Scenario = mf$Scenario,
  SpeciesPool = mf$SpeciesPool,
  ForestID = mf$ForestID
)

mf <- model.frame(mem_iqr)            # model frame used to fit mem_iqr
df <- data.frame(
  resid        = resid(mem_iqr, type = "pearson"),
  fitted       = fitted(mem_iqr),
  Scenario     = mf$Scenario,
  SpeciesPool  = mf$SpeciesPool,
  ForestID     = mf$ForestID
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

pdf(file.path(DirectoryPlots, "pos_res_grouped_v2.pdf"), width = 7, height = 7)
print(p2 / p1)
dev.off()

# -- Normality of randm effects

# Extract random effects
species_pool <- ranef(mem_pos)$SpeciesPool[,1]
forest_id <- ranef(mem_pos)$ForestID[,1]

# Plot histograms
pdf(file.path(DirectoryPlots, "diag_pos_rf_sp_normal_glmer.pdf"))
hist(species_pool)
dev.off()

pdf(file.path(DirectoryPlots, "diag_pos_rf_forest_normal_glmer.pdf"))
hist(forest_id)
dev.off()

# Test for normality
shapiro.test(species_pool)
shapiro.test(forest_id)

# Variance partitioning

icc(mem_pos) # -> random effects explain <1% of variance

# Variance components
sp_var <- 0.0003057
frst_var <- 0.0009625
total_var <- sp_var + frst_var + sigma(mem_pos)^2
sp_vc <- 100*(sp_var / total_var)
frst_vc <- 100*(frst_var / total_var)

print(paste0("Species pool variance component: ", round(sp_vc, 4)))
print(paste0("Forest ID variance component: ", round(frst_vc, 4)))


r2 <- performance::r2_nakagawa(mem_pos)
icc <- performance::icc(mem_pos)
icc_unadj <- icc$ICC_unadjusted * 100

# Get variances
vc <- VarCorr(mem_pos)

# Get variances
var_species_intercept <- vc$SpeciesPool["(Intercept)", "(Intercept)"]
var_forest             <- vc$ForestID["(Intercept)", "(Intercept)"]

# Total random-effect variance
total_var <- var_species_intercept + var_forest

# Compute proportions
prop_species_intercept <- (var_species_intercept / total_var) * icc_unadj
prop_forest            <- (var_forest / total_var) * icc_unadj

cat(
  paste0(
    "Conditional R2: ", round(r2$R2_conditional, 3), "\n",
    "Marginal R2: ", round(r2$R2_marginal, 3), "\n",
    "ICC unadjusted: ", round((icc$ICC_unadjusted / r2$R2_conditional) * 100, 2), "%\n",
    "Sp: ", round(var_species_intercept, 2), "(", round(prop_species_intercept, 2), "%)\n",
    "Forest: ", round(var_forest, 2), "(", round(prop_forest, 2), "%)"
  )
)
