library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(permute)

numSpeciesPools <- seq(1, 10)
replicatePerSpeciesPool <- 1
forestIDs <- c(0, 1, 2)
timeStepStart <- 80
timeStepEnd <- 199

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

  DirectoryModelResults <- paste0("/Users/johanna/Uni/masterarbeit/data/modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/a5/forest", forestId, "/")
  DirectoryPlots <- paste0("../../../figs/a5_plots_test/climdata_era5_cmip6_1981-2100_ssp245/forest", forestId, "/")

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

write_csv(
  speciesHeight,
  "/Users/johanna/Uni/masterarbeit/data/modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/a5/species_heights.csv",
)

overallStats <- speciesHeight %>%
  group_by(TimeStep, Replicate, SpeciesPool, ForestID) %>%
  summarize(
    Richness = length(unique(SpeciesID)),
    Abundance = n()
  )

# ---- 1. approach

set.seed(42)

n_forests <- 3
n_species_pools <- 10
n_time_steps <- 120

for (i in 1:10) {
  results_list <- list()
  shuffled_reps <- shuffle(1:30)

  for (draw_size in seq(3, 30, 1)) {
    samples <- shuffled_reps[1:draw_size]
    # convert to forest + species
    idx <- data.frame(arrayInd(samples, .dim = c(n_forests, n_species_pools)))
    colnames(idx) <- c("ForestID", "SpeciesPool")

    for (t in 80:199) {   # loop over time steps
      # join with overallStats
      dat <- idx %>%
        mutate(
          ForestID = ForestID - 1,
          SpeciesPool = SpeciesPool,
          TimeStep = t,
          Replicate = 1
        ) %>%
        left_join(overallStats,
                  by = c("ForestID", "SpeciesPool", "TimeStep", "Replicate"))

      # compute CVs
      richnessCV <- sd(dat$Richness) / mean(dat$Richness)
      abundanceCV <- sd(dat$Abundance) / mean(dat$Abundance)

      results_list[[length(results_list) + 1]] <- data.frame(
        TimeStep = t,
        SampleSize = draw_size,
        RichnessCV = richnessCV,
        AbundanceCV = abundanceCV
      )
    }
  }

  CVs <- bind_rows(results_list) %>% tibble()

  pRichness <- ggplot(CVs, aes(x = TimeStep, y = RichnessCV, group = SampleSize, color = factor(SampleSize))) +
  geom_line(alpha = 0.3) +
  labs(
    x = "Time Step",
    y = "Richness CV"
  ) +
  theme_minimal()

pAbundance <- ggplot(CVs, aes(x = TimeStep, y = AbundanceCV, group = SampleSize, color = factor(SampleSize))) +
  geom_line(alpha = 0.3) +
  labs(
    x = "Time Step",
    y = "Abundance CV"
  ) +
  theme_minimal()

  pdf(paste0("../../../figs/a5_plots_test/", i, "CV_bootstrap_richness_abundance_incremental_cc_ts_v3.pdf"),
      width = 12, height = 6)
  print(pRichness + pAbundance)
  dev.off()
}

# ---- 2. approach

set.seed(42)

community_last10ts <-  overallStats %>%
    filter(TimeStep <= 199 & TimeStep >= 189) %>%
    group_by(Replicate, SpeciesPool, ForestID) %>%
    summarise(
      Richness = mean(Richness, na.rm = TRUE),
      Abundance = mean(Abundance, na.rm = TRUE)
    )

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

# plots
pRichness <- ggplot(CVsLast10, aes(x = SampleSize, y = RichnessCV, group = Rep)) +
  geom_line(alpha = 0.3, color = "steelblue") +
  labs(
    x = "Number of Draws",
    y = "Richness CV"
  ) +
  theme_minimal()

pAbundance <- ggplot(CVsLast10, aes(x = SampleSize, y = AbundanceCV, group = Rep)) +
  geom_line(alpha = 0.3, color = "darkred") +
  labs(
    x = "Number of Draws",
    y = "Abundance CV"
  ) +
  theme_minimal()

pdf("../../../figs/a5_plots_test/CV_bootstrap_richness_abundance_no_ts_cc_v3.pdf", width = 12, height = 6)
print(pRichness + pAbundance)
dev.off()


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

pdf("../../../figs/a5_plots_test/CV_bootstrap_richness_abundance_no_ts_cc_mean_sd.pdf",
    width = 12, height = 6)
print(pRichnessMean + pAbundanceMean)
dev.off()