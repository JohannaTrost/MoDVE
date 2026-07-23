# -----
# Generate dataset of vertical species distribution combining all MC gradient scenarios

library(dplyr)
library(tidyr)
library(vegan)
library(readr)

DirectoryModelResults <- file.path("../modve_data/modve_output/pirineus/scenarios")
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

speciesVertical <- NULL

for (scenario in names(senarios)) {
  cat(scenario, "\n")
  # Loop over species pools
  for (numSpeciesPool in numSpeciesPools) {
    cat(numSpeciesPool, "\n")
    # Collecting data for each time step
    for (timeStep in seq(timeStepStart, timeStepEnd)) {

      modelResultsPath <- file.path(
          DirectoryModelResults, scenario, "communities",
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

      # - Species verticla metrics
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
      if (is.null(speciesVertical)) {
        speciesVertical <- speciesVerticalCurr
      } else {
        # combine new and old
        speciesVertical <- speciesVerticalCurr %>%
          rbind(., speciesVertical)
      }
    }
  }
}

savePathSpec <- file.path(
          DirectoryModelResults,
          "SpeciesVertical_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_mcGradients.csv"
)
write_csv(speciesVertical, savePathSpec)