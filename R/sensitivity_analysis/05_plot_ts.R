library(dplyr)
library(readr)
#library(ggsignif)
library(ggpubr)
library(ggplot2)

###############################################################################################
#                                        Diversity data                                       #
###############################################################################################

DirectoryModelResults <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")

divFile1 <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxRichness_mcGradients.csv"
divFile2 <- "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxShannonDiv_mcGradients.csv"

maxRichness <- read_csv(file.path(DirectoryModelResults, divFile1), show_col_types = FALSE)
maxShannon <- read_csv(file.path(DirectoryModelResults, divFile2), show_col_types = FALSE)

############# Richness peak position #############

# Avg over space
avg_richness_pos <- maxRichness %>%
  group_by(scenario, speciesPool, timeStep) %>%
  summarise(meanMaxZ = mean(maxZ), .groups = "drop")

# Summarize across species pools for each scenario and timestep
avg_summary <- avg_richness_pos %>%
  group_by(scenario, timeStep) %>%
  summarise(
    mean_richness = mean(meanMaxZ, na.rm = TRUE),
    sd_richness   = sd(meanMaxZ, na.rm = TRUE),
    .groups = "drop"
  )

# - Time series Plot
pdf("../../figs/sensitivity_analysis/maxRichnessZ_mcGradients_ts.pdf")
ggplot(avg_summary, aes(x = timeStep, y = mean_richness,
                        group = scenario, linetype = factor(scenario))) +
  geom_ribbon(aes(ymin = mean_richness - sd_richness,
                  ymax = mean_richness + sd_richness,
                  fill = factor(scenario)), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  labs(
    x = "Time step",
    y = "Avg. vertical richness peak position (m)",
    linetype = "Scenario",
    fill = "Scenario"
  ) +
  theme_minimal(base_size = 14)
dev.off()

# - Box plot
maxTs <- max(maxRichness$timeStep)

pdf("../../figs/sensitivity_analysis/maxRichnessZ_mcGradients_boxplot_v3.pdf")
maxRichness %>%
  filter(timeStep >= (maxTs - 30)) %>%
  ggplot(aes(x = factor(scenario), y = maxZ)) +
    geom_boxplot(notch = TRUE, aes(fill = factor(scenario), alpha = 0.2)) +
    labs(
      x = "Steepness factor of MC gradient",
      y = "Vertical richness peak position (m)"
    ) +
    theme_minimal() +
    # Pairwise comparisons with multiple testing correction
    stat_compare_means(
      comparisons = list(c("0", "0.5"), c("0", "1"), c("0.5", "1"), c("1.5", "1")), # adjust as needed
      method = "wilcox.test",
      p.adjust.method = "bonferroni", # or "fdr", "holm"
      label = "p.signif"
    ) +
    stat_compare_means( # global test across groups
      method = "anova",
      label.y = max(maxRichness$maxZ) * 1.1
    )
dev.off()

# - Effect size
eff_sizes <- lapply(pairs, function(p) {
  df_sub <- maxRichness %>%
    filter(timeStep >= (maxTs - 30),
           scenario %in% p) %>%
    as.data.frame()

  # Cohen's d using formula
  d_res <- cohens_d(maxZ ~ factor(scenario), data = df_sub)

  data.frame(
    group1 = p[1],
    group2 = p[2],
    cohen_d = d_res$Cohens_d,
    conf_low = d_res$CI_low,
    conf_high = d_res$CI_high
  )
})

eff_sizes <- bind_rows(eff_sizes)
eff_sizes

# TODO test normality

# TODO do this for richnes itself

# TODO do this for species position

# TODO do this for species range