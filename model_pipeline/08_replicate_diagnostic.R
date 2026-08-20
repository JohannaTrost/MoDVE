# -----
# Evaluate whether the number of replicates is sufficient with bootstrapping approach:
# 1. Draw an increasing number of replicates (1 to 30)
# 2. Compute the coefficient of variation (CV) across richness and abundance of these replicates
# 3. Repeat 100 times
# 4. Plot the number of draws against the CV

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(permute)

# Directories
base_dir <- file.path("../modve_output/regua")
figs_dir <- file.path("../modve_figs/climdata_era5_cmip6_1981-2100_ssp245")

# Parameters
numSpeciesPools <- seq(1, 10)
replicatePerSpeciesPool <- 1
forestIDs <- c(0, 1, 2)
timeStepStart <- 80
timeStepEnd <- 199

set.seed(42) # Seed for reproducability

speciesHeight <- tibble(
  SpeciesPool = integer(),
  Replicate = integer(),
  TimeStep = integer(),
  SpeciesID = integer(),
  height = numeric(),
  ForestID = integer()
)

for (forestId in forestIDs) {

  print(paste0("Processing Forest: ", forestId))

  DirectoryModelResults <- file.path(
    base_dir, "climdata_era5_cmip6_1981-2100_ssp245", "communities", paste0("forest", forestId)
  )

  for (rep in seq_len(replicatePerSpeciesPool)) {
    for (numSpeciesPool in numSpeciesPools) {

      print(paste0("Species Pool: ", numSpeciesPool, " Replicate: ", rep))

      # --- Collecting data for each time step
      for (timeStep in int_seq(timeStepStart, timeStepEnd)) {

        modelResultsPath <- paste0(DirectoryModelResults,
                                   "ID_SpeciesP_", numSpeciesPool, "_Rep_", rep,
                                   "/IndividualMatrixTimeStep", timeStep, ".csv")
        if (!file.exists(modelResultsPath)) {
          message("File does not exist: ", modelResultsPath)
          next
        }

        # Load results
        res <- read_csv(modelResultsPath, show_col_types = FALSE)
        res <- res[res$Status == 1,]  # We are only interested in individuals that survived

        if (nrow(res) == 0) {
          message("No individuals survived at time step ", timeStep, " for Species Pool ", numSpeciesPool, " and Replicate ", rep)
          next
        }

        # Collect species height data
        res["height"] <- res$Z - 0.5
        speciesHeight <- speciesHeight %>%
          add_row(SpeciesPool = numSpeciesPool,
                  Replicate = rep,
                  TimeStep = timeStep,
                  SpeciesID = res$SpeciesID,
                  height = res$height,
                  ForestID = forestId)
      }
    }
  }
}

# Get richness and abundance
overallStats <- speciesHeight %>%
  group_by(TimeStep, Replicate, SpeciesPool, ForestID) %>%
  summarize(
    Richness = length(unique(SpeciesID)),
    Abundance = n()
  )

# Get status across last 10 time steps
community_last10ts <-  overallStats %>%
    filter(TimeStep <= 199 & TimeStep >= 189) %>%
    group_by(Replicate, SpeciesPool, ForestID) %>%
    summarise(
      Richness = mean(Richness, na.rm = TRUE),
      Abundance = mean(Abundance, na.rm = TRUE)
    )

# --- 100 iterations of Bootstrap sampling

results_list <- list()

for (i in 1:100) {

  shuffled_reps <- shuffle(1:30)

  for (draw_size in seq(3, 30, 1)) {
    samples <- shuffled_reps[1:draw_size]
    # convert to forest + species
    idx <- data.frame(arrayInd(samples, .dim = c(n_forests, n_species_pools)))
    colnames(idx) <- c("ForestID", "SpeciesPool")

    dat <- idx %>%
      mutate(
        ForestID = ForestID - 1,
        SpeciesPool = SpeciesPool,
        Replicate = 1
      ) %>%
      left_join(community_last10ts,
                by = c("ForestID", "SpeciesPool", "Replicate"))

    # compute CVs
    richnessCV <- sd(dat$Richness) / mean(dat$Richness)
    abundanceCV <- sd(dat$Abundance) / mean(dat$Abundance)

    results_list[[length(results_list) + 1]] <- data.frame(
      SampleSize = draw_size,
      RichnessCV = richnessCV,
      AbundanceCV = abundanceCV,
      Rep = i
    )
  }
}

CVsLast10 <- bind_rows(results_list) %>% tibble()

# --- Plot the results

pRichness <- ggplot(CVsLast10, aes(x = SampleSize, y = RichnessCV, group = Rep)) +
  geom_line(alpha = 0.3, color = "steelblue") +
  labs(
    x = "Number of Draws",
    y = "Richness CV"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 25),
    axis.text = element_text(size = 20))

pAbundance <- ggplot(CVsLast10, aes(x = SampleSize, y = AbundanceCV, group = Rep)) +
  geom_line(alpha = 0.3, color = "darkred") +
  labs(
    x = "Number of Draws",
    y = "Abundance CV"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 25),
    axis.text = element_text(size = 20))

pdf(file.path(figs_dir, "CV_bootstrap_richness_abundance_no_ts_cc_v4.pdf"), width = 12, height = 6)
print(pRichness + pAbundance)
dev.off()

# --- Plot mean and SD across iterations

# Plot average and +-SD as ribbon over time steps but in principles samplot as above
pRichnessMean <- ggplot(CVsLast10, aes(x = SampleSize, y = RichnessCV)) +
  stat_summary(fun = mean, geom = "line", color = "steelblue") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  labs(
    x = "Number of Draws",
    y = "Richness CV"
  ) +
  theme_minimal()

pAbundanceMean <- ggplot(CVsLast10, aes(x = SampleSize, y = AbundanceCV)) +
  stat_summary(fun = mean, geom = "line", color = "coral") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "ribbon", alpha = 0.2, fill = "coral") +
  labs(
    x = "Number of Draws",
    y = "Abundance CV"
  ) +
  theme_minimal()

pdf(file.path(figs_dir, "CV_bootstrap_richness_abundance_no_ts_cc_mean_sd.pdf"), width = 12, height = 6)
print(pRichnessMean + pAbundanceMean)
dev.off()