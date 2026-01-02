library(dplyr)
library(tidyr)
library(vegan)
library(readr)


DirectoryModelResults <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")
numSpeciesPools <- seq(1, 10)
replicatePerSpeciesPool <- 1
timeStepStart <- 80
timeStepEnd <- 197
stepSize <- 5

senarios <- list(
  "climdata_era5_cmip6_1906-2024_ssp245_119ts_no_mc_grad" = 0,
  "climdata_era5_cmip6_1906-2024_ssp245_119ts_0.5_mc_grad" = 0.5,
  "climdata_era5_cmip6_1906-2024_ssp245_119ts" = 1,
  "climdata_era5_cmip6_1906-2024_ssp245_119ts_1.5_mc_grad" = 1.5
)

diversity <- NULL

for (scenario in names(senarios)) {
  cat(scenario, "\n")
  # Loop over species pools
  for (numSpeciesPool in numSpeciesPools) {
    cat(numSpeciesPool, "\n")
    # Collecting data for each time step
    for (timeStep in seq(timeStepStart, timeStepEnd)) {

      modelResultsPath <- file.path(
          DirectoryModelResults, scenario, "a5", "spec200_rand",
          paste0("ID_SpeciesP_", numSpeciesPool, "_Rep_1/IndividualMatrixTimeStep",
                 timeStep, ".csv"))
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

      # - 1. Diveristy and richness
      # Compute species abundances per voxel
      voxel_abund <- res %>%
        group_by(X, Y, Z, SpeciesID) %>%
        summarise(n = n(), .groups = "drop") %>%   # counts individuals per species per voxel
        pivot_wider(
          names_from = SpeciesID,
          values_from = n,
          values_fill = 0
        )

      # Extract voxel coordinates separately
      coords <- voxel_abund %>% select(X, Y, Z)
      abund_mat <- voxel_abund %>% select(-X, -Y, -Z)

      # Compute richness and Shannon diversity per voxel
      richness <- specnumber(abund_mat)   # number of species per voxel
      shannon  <- diversity(abund_mat, index = "shannon")  # Shannon entropy

      # Merge res into df
      curr_div <- coords %>%
          mutate(
            richness = richness,
            shannon = shannon,
            timeStep = timeStep,
            speciesPool = numSpeciesPool,
            scenario = as.numeric(senarios[scenario])
          )

      # - 2. Species verticla metrics

      speciesVerticalCurr <- res %>%
         mutate(
           Zcenter = Z - 0.5   # voxel center height
         ) %>%
         group_by(SpeciesID) %>%
         summarise(
           timeStep = timeStep,
           speciesPool = numSpeciesPool,
           scenario = as.numeric(senarios[scenario]),
           meanHeight = mean(Zcenter, na.rm = TRUE),
           varHeight  = var(Zcenter, na.rm = TRUE),
           iqrHeight  = IQR(Zcenter, na.rm = TRUE),
           n_individuals = n(),         # number of individuals per species (for context)
           .groups = "drop"
         )

            # Combine results
      if (is.null(diversity)) {
        diversity <- curr_div
        speciesVertical <- speciesVerticalCurr
      } else {
        # combine new and old
        diversity <- curr_div %>%
          rbind(., diversity)
        speciesVertical <- speciesVerticalCurr %>%
          rbind(., speciesVertical)
      }
    }
  }
}

# --- Complete misisng values with 0s for richness and shannon

# Define the full set of combinations
full_combinations <- expand.grid(
  X = unique(diversity$X),
  Y = unique(diversity$Y),
  Z = unique(diversity$Z),
  speciesPool = numSpeciesPools,
  scenario = unlist(senarios),
  timeStep = seq(timeStepStart, timeStepEnd, by = stepSize)
)

# Join with existing data
diversity_complete <- full_combinations %>%
  left_join(diversity, by = c("X", "Y", "Z", "speciesPool", "scenario", "timeStep")) %>%
  # Replace NA richness and shannon with 0
  mutate(
    richness = ifelse(is.na(richness), 0, richness),
    shannon = ifelse(is.na(shannon), 0, shannon)
  ) %>%
  # Remove duplicates just in case
  distinct(X, Y, Z, speciesPool, scenario, timeStep, .keep_all = TRUE)

# Compute max richness and max shannon diversity per (X,Y) location, species pool, time step, scenario
maxRichness <- diversity_complete %>%
  group_by(X, Y, speciesPool, timeStep, scenario) %>%
  slice_max(order_by = richness, with_ties = FALSE) %>%
  transmute(
    X, Y, speciesPool, timeStep, scenario,
    maxZ = Z,
    maxRichness = richness
  ) %>%
  ungroup()

# Not by cell
richness_by_Z <- diversity_complete %>%
  group_by(speciesPool, timeStep, scenario, Z) %>%
  summarise(
    richness_sum = sum(richness, na.rm = TRUE),
    .groups = "drop"
  )

maxRichnessAgg <- richness_by_Z %>%
  group_by(speciesPool, timeStep, scenario) %>%
  slice_max(order_by = richness_sum, with_ties = FALSE) %>%
  rename(
    maxRichness = richness_sum,
    maxZ = Z
  ) %>%
  ungroup()

maxShannon <- diversity_complete %>%
  group_by(X, Y, speciesPool, timeStep, scenario) %>%
  slice_max(order_by = shannon, with_ties = FALSE) %>%
  transmute(
    X, Y, speciesPool, timeStep, scenario,
    maxZ = Z,
    maxDiversity = shannon
  ) %>%
  ungroup()

# Save results
savePathDiv1 <- file.path(
          DirectoryModelResults,
          "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxRichness_mcGradients_full.csv")
write_csv(maxRichness, savePathDiv1)
savePathDiv1 <- file.path(
          DirectoryModelResults,
          "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxRichnessAgg_mcGradients_full.csv")
write_csv(maxRichnessAgg, savePathDiv1)
savePathDiv1 <- file.path(
          DirectoryModelResults,
          "a5_Rep_1_climdata_era5_1906-2024_119ts_verticalRichness_mcGradients_full.csv")
write_csv(richness_by_Z, savePathDiv1)
savePathDiv2 <- file.path(
          DirectoryModelResults,
          "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_maxShannonDiv_mcGradients_full.csv")
write_csv(maxShannon, savePathDiv2)
savePathSpec <- file.path(
          DirectoryModelResults,
          "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_SpeciesVertical_mcGradients.csv"
)
write_csv(speciesVertical, savePathSpec)