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

mem_pos <- glmmTMB(
  Range ~ Scenario * Year_c + (Scenario | SpeciesPool) + (1 | ForestID),
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

# QQ plot
pdf(file.path(DirectoryPlots, "diag_range_res_qq_dharma_glmmtmb2_v5.pdf"))

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
