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

# #################################################################################################
#                               Why do some species shift upward?                                 #
# #################################################################################################

# Compute position shift
shift <- species_distr_stats %>%
  dplyr::select(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year, Position) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = Position
  ) %>%
  mutate(
    diff = `CC` - `No CC`
  ) %>%
  filter(!is.na(diff))

summary_shift <- shift %>%
  group_by(ForestID, SpeciesPool, SpeciesID) %>%
  arrange(Year, .by_group = TRUE) %>%  # ensure data is sorted by Year
  summarise(
    initial_diff = diff[Year == 2000],
    diff = mean(tail(diff[!is.na(diff)], 10), na.rm = TRUE),
    last_year_alive = max(Year[!is.na(diff)]),
    .groups = "drop"
  )

# Get species niches
niches <- NULL
for (sp in 1:10) {
  sp_niches <- read_csv(file.path(base_dir, "a2_1", paste0("SpeciesPool", sp, ".csv")),
                        show_col_types = FALSE) %>%
    dplyr::select("SpeciesID", "OptimumLight", "OptimumTemp", "OptimumHum",
           "MinLight", "MaxLight", "MinTemp", "MaxTemp", "MinHum", "MaxHum") %>%
    mutate(SpeciesPool = sp)

  if (is.null(niches)) {
    niches <- sp_niches
  } else {
    niches <- rbind(niches, sp_niches)
  }
}

# Compute ranges
niches$RangeLight <- niches$MaxLight - niches$MinLight
niches$RangeTemp <- niches$MaxTemp - niches$MinTemp
niches$RangeHum <- niches$MaxHum - niches$MinHum

# Combine with shift data
shift_niches <- summary_shift %>%
  left_join(niches, by = c("SpeciesPool", "SpeciesID"))

# Compute species distances
cols_to_use <- c("OptimumLight", "OptimumTemp", "OptimumHum",
                 "MinLight", "MaxLight", "MinTemp", "MaxTemp", "MinHum", "MaxHum",
                 "RangeLight", "RangeTemp", "RangeHum")

shift_niches_dist <- shift_niches %>%
  group_by(SpeciesPool, ForestID) %>%
  mutate(across(all_of(cols_to_use),
                ~ {
                    col_values <- .
                    sapply(seq_along(col_values), function(i) {
                      mean(abs(col_values[i] - col_values[-i]))
                    })
                  },
                .names = "MeanDist_{.col}")) %>%
  ungroup()

# 1. Visualize with PCA

metric <- c("Min", "Max")

# Select all MeanDist columns for PCA
dist_cols <- c()
for (m in metric) {
  dist_cols <- c(dist_cols,
                 grep(paste0("^MeanDist_", m), names(shift_niches_dist), value = TRUE))
}
pca_data <- shift_niches_dist[, dist_cols]

# Remove any rows with missing values
pca_data_complete <- na.omit(pca_data)
complete_indices <- complete.cases(shift_niches_dist[, dist_cols])
data_for_plot <- shift_niches_dist[complete_indices, ]

# Perform PCA
pca_result <- prcomp(pca_data_complete, scale. = TRUE, center = TRUE)

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Print variance explained for PC1 and PC2
cat("Variance Explained:\n")
cat(sprintf("PC1: %.2f%%\n", var_explained[1]))
cat(sprintf("PC2: %.2f%%\n", var_explained[2]))
cat(sprintf("Total (PC1 + PC2): %.2f%%\n", sum(var_explained)))

# Create dataframe for plotting
plot_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  SpeciesID = data_for_plot$SpeciesID,
  diff = data_for_plot$diff,
  positive_diff = data_for_plot$diff > 0
)

# Create the plot
p <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = positive_diff), alpha = 0.6, size = 2) +
  geom_text(data = subset(plot_df, positive_diff),
            aes(label = SpeciesID),
            size = 3,
            hjust = -0.1,
            vjust = 0.5,
            check_overlap = TRUE) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"),
                     labels = c("FALSE" = "diff ≤ 0", "TRUE" = "diff > 0"),
                     name = "Difference") +
  labs(
    title = "PCA of Niche Distance Metrics",
    subtitle = sprintf("PC1: %.2f%% | PC2: %.2f%% variance explained",
                      var_explained[1], var_explained[2]),
    x = sprintf("PC1 (%.2f%%)", var_explained[1]),
    y = sprintf("PC2 (%.2f%%)", var_explained[2])
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom"
  )

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/position_shift")
dir.create(DirectoryPlots)

pdf(file.path(DirectoryPlots, paste0("pc_distances_", paste(metric, collapse = "_"), ".pdf")))
print(p)
dev.off()

# --- Try MEM

pdf(file.path(DirectoryPlots, "diag_shift_hist.pdf"))
hist(shift_niches_dist$diff)
dev.off()

shapiro.test(shift_niches_dist$diff)

# Scale predictors
shift_niches_dist <- shift_niches_dist %>%
  mutate(across(starts_with("MeanDist_"),
                ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),
                .names = "Scaled_{.col}"))

lm_shift <- lm(
   diff ~ MeanDist_RangeTemp *
     MeanDist_RangeHum *
     MeanDist_RangeLight,
  data = shift_niches_dist)

summary(lm_shift)