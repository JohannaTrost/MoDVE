
run_modve_sim <- function(sim_params,
                          SpeciesPool,
                          Microhabitat,
                          InitDist,
                          path_to_ind_output,
                          path_to_sp_output,
                          path_to_comm_output) {

  # Check parameters and read input if necessary
  exptd_params <- c(
    "InitialTimeStep", "timeSteps", "StopCriterionHa","MicrohabitatType",
    "Imax", "CompetitionMethod", "MortalityMethod", "MortRateMass",
    "MortRateMassScaling", "MortRateRandom", "SurfaceBiomassScaling",
    "SlopeRecruitment", "InterceptRecruitment"
  )
  if (!all(is.element(exptd_params, sim_params))) {
    stop(cat("sim_params must contain at least the following elements: ",
             exptd_params, "\n"))
  }
  timeSteps <- sim_params$timeSteps

  if (!is.data.frame(SpeciesPool)) {
    # Then it must be a path to the input
    if (!grepl("*.csv$", SpeciesPool)) {
      stop("SpeciesPool should be a data frame or a valid path to a .csv file.")
    }
    if (!file.exists(SpeciesPool)) {
      stop(paste0(SpeciesPool, " doesn't exist.\n"))
    } # error
    else SpeciesPool <- read.csv(SpeciesPool, sep = ",", header = TRUE)
  }
  NumberOfSpecies <- nrow(SpeciesPool)  # number of species per 25X25m plot

  if (!is.data.frame(InitDist)) {
    # Then it must be a path to the input
    if (!grepl("*.csv$", InitDist)) {
      stop("InitDist should be a data frame or a valid path to a .csv file.")
    }
    if (!file.exists(InitDist)) {
      stop(paste0(InitDist, " doesn't exist.\n"))
    }
    else InitDist <- read.csv(InitDist, sep = ",", header = TRUE)
  }
  E <- InitDist
  # Add columns to E for additional info
  E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0
  max_id <- nrow(E)  # to trace individual IDs

  isHabitatDynamic <- sim_params$MicrohabitatType == 1
  if (!is.data.frame(Microhabitat)) {
    # Then it must be a path or vector of paths
    for (i in seq_along(Microhabitat[i])) {
      if (!grepl("*.csv$", Microhabitat[i])) {
        stop("Microhabitat should be a data frame or a valid path to a .csv file.")
      }
      if (!file.exists(Microhabitat[i])) {
        stop(paste0(Microhabitat[i], " doesn't exist.\n"))
      }
    }

    if (isHabitatDynamic) {
      if (!length(Microhabitat) == timeSteps) {
        stop("For dynamic habitats, Microhabitat should have one element for each time step.\n")
      } else {
        # Stash paths for later time steps
        microhab_files <- Microhabitat
      }
    }
    # If all checks ok, read the first one
    Microhabitat <- read.csv(Microhabitat[1], sep = ",", header = TRUE)
  }

  #  Convert relative light values to absolute ?mol*m-2*s-1
  Microhabitat[, , , 3] <- Microhabitat[, , , 3] * sim_params$Imax
  dimX <- dim(Microhabitat)[1]
  dimY <- dim(Microhabitat)[2]
  dimZ <- dim(Microhabitat)[3]

  # Prepare output
  dir_output <- dirname(path_to_ind_output)
  if (!dir.exists(dirname(dir_output))) {
    stop(paste0("Directory ", dir_output, " does not exist."))
  }
  if (!grepl("*.csv$", path_to_ind_output)) {
    stop("path_to_ind_output must be a csv file")
  }

  dir_output <- dirname(path_to_sp_output)
  if (!dir.exists(dirname(dir_output))) {
    stop(paste0("Directory ", dir_output, " does not exist."))
  }
  if (!grepl("*.csv$", path_to_sp_output)) {
    stop("path_to_sp_output must be a csv file")
  }


  sp_output_headers <- species_output_names()
  {
    # Column indices
    col_sp_t <- 1; col_sp_id <- 2; col_nb_inds_begin <- 3; col_nb_inds_end <- 4;
    col_nb_mature_inds <- 5; col_nb_rec <- 6; col_nb_rec <- 7;
    col_nb_dead_branch <- 8; col_nb_dead_light <- 9; col_nb_dead_comp <- 10;
    col_nb_dead_base <- 11; col_growth_rate <- 12; col_growth_log <- 13
    col_birth <- 14; col_death <- 15; col_size <- 16; ColSAverageAge <- 17
    ColSMinLight <- 18; ColSMaxLight <- 19; ColSMeanLight <- 20;
    ColSMinHeight <- 21; ColSMaxHeight <- 22; ColSMeanHeight <- 23
    nb_cols_sp_output <- ColSMeanHeight
  }
  # Initialize Matrix where species parameters are saved
  sp_output <- array(
    rep(0, (timeSteps * NumberOfSpecies) * nb_cols_sp_output),
    dim = c(timeSteps * NumberOfSpecies, nb_cols_sp_output)
  )

  dir_output <- dirname(path_to_comm_output)
  if (!dir.exists(dirname(dir_output))) {
    stop(paste0("Directory ", dir_output, " does not exist."))
  }
  if (!grepl("*.csv$", path_to_comm_output)) {
    stop("path_to_comm_output must be a csv file")
  }
  comm_output_headers <- comm_output_names()
  comm_output <- data.frame(matrix(
    0.0, nrow = timeSteps, ncol = length(comm_output_headers)
  ))
  colnames(comm_output) <- comm_output_headers

  # Stop criterion: if stop density is exceeded, the simulation ends
  # to cap memory usage
  StopNbInds <- sim_params$StopCriterionHa * dimX * dimY / 10000

  # Calculate the probability to disperse in surrounding voxels
  centralPoint <- find_central_point(c(dimX, dimY, dimZ))
  prob_disp_matrix <- calc_prob_disp_matrix(
    centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool
  )

  # Generation loop
  for (t in seq_len(timeSteps)) {

    # Check if the stop criterion is met
    nbIndsAlive <- length(which(E$Status == 1))
    if (nbIndsAlive > StopNbInds) {
      writeLines(paste0(
        "Time ", t, ": population has exceeded max threshold of ",StopNbInds,
        ". Ending simulation."
        ))
      break
    }

    # Update microhabitat if applicable
    if (isHabitatDynamic && t > 1) {
      Microhabitat <- readRDS(microhab_files[t])
      Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax
    }

    # Update how many species are alive at beginning of generation
    InitialNumberSpecies <- length(unique(E$SpeciesID[E$Status == 1]))
    nbIndsBeforeDispTotal <- length(which(E$Status == 1))

    # Dispersal
    disp_items <- resolveReproDispersal(
      E, Microhabitat, SurfaceBiomassScaling,
      centralPoint, InterceptRecruitment, SlopeRecruitment,
      prob_disp_matrix,  SpeciesPool, max_id
    )

    # Unwrap dispersal output
    nbIndsBeforeDisp <- disp_items$nbIndsBeforeDisp
    recruitment_df <- disp_items$recruitment_df
    E <- disp_items$E
    max_id <- disp_items$max_id

    # Store potential normalized number of recruits in sp_output
    for (i in seq_len(nrow(recruitment_df))) {
      sp_idx <- recruitment_df$index[i]
      row_idx <- (sp_idx - 1) * timeSteps + t
      sp_output[row_idx, col_nb_rec] <-
        recruitment_df$nb_potential_recruits[i]
    }
    NumberRecruits <- length(which(E$Status == 1)) - nbIndsBeforeDispTotal
    nbRecruitsPerSpecies <- disp_items$recruitment_df$nb_recruits

    # TODO: Unclear what this line in the Matlab script is supposed to do.
    # From what I understand, the first column in E ("SpeciesID") takes non-zero values
    # only, so I think that E(:,1)==0 will always be empty.
    # E(E(:,1)==0,:)=[]; %in rare case, some individuals with only zeros are creates, which is wrong. This is to prevent the script to stop.

    # Growth
    E <- resolve_growth(E, Microhabitat, SurfaceBiomassScaling)

    # Mortality (except from competition)
    E <- resolve_mortality(E, Microhabitat)

    # Mortality due to competition for space
    E <- resolve_competition(E, Microhabitat, CompetitionMethod)

    # Age increment
    E$Age <- E$Age + 1

    # Species-level output
    for (nb_sp in seq_len(NumberOfSpecies)) {

      nb_alive <- sum(E$Status == 1 & is_sp, na.rm = TRUE)
      nb_dead_comp <- sum(E$Status == 2 & is_sp, na.rm = TRUE)
      nb_dead_branch <- sum(E$Status == 3 & is_sp, na.rm = TRUE)
      nb_dead_light <- sum(E$Status == 4 & is_sp, na.rm = TRUE)
      nb_dead_base <- sum(E$Status == 5 & is_sp, na.rm = TRUE)
      nb_inds_begin <- nbIndsBeforeDisp[nb_sp]

      row_nb <- (nb_sp - 1) * timeSteps + t
      is_sp <- E$SpeciesID == nb_sp

      sp_output[row_nb, col_sp_t] <- InitialTimeStep + t - 1
      sp_output[row_nb, col_sp_id] <- nb_sp
      sp_output[row_nb, col_nb_inds_begin] <- nb_inds_begin
      sp_output[row_nb, col_nb_inds_end] <- nb_alive
      sp_output[row_nb, col_nb_mature_inds] <- sum(
        E$Status == 1 & is_sp & E$Mass >= E$MassAtMaturity,
        na.rm = TRUE
      )
      sp_output[row_nb, col_nb_rec] <- nbRecruitsPerSpecies[nb_sp]
      sp_output[row_nb, col_nb_dead_branch] <- nb_dead_branch
      sp_output[row_nb, col_nb_dead_light] <- nb_dead_light
      sp_output[row_nb, col_nb_dead_comp] <- nb_dead_comp
      sp_output[row_nb, col_nb_dead_base] <- nb_dead_base

      if (nb_alive > 0 &&  nb_inds_begin > 0) {
        sp_output[row_nb, col_growth_rate] <- sp_output[row_nb, col_nb_inds_end] /
          sp_output[row_nb, col_nb_inds_begin]
        sp_output[row_nb, col_growth_log] <- log(sp_output[row_nb, col_growth_rate])
        sp_output[row_nb, col_birth] <- nbRecruitsPerSpecies[nb_sp] /
          nb_inds_begin
        sp_output[row_nb, col_death] <- (
          nb_dead_branch + nb_dead_light + nb_dead_comp + nb_dead_base
        ) / nb_inds_begin
        sp_output[row_nb, col_size] <- mean(E$Mass[is_sp])
        sp_output[row_nb, ColSAverageAge] <- mean(E$Age[is_sp])
        sp_output[row_nb, ColSMinLight] <- min(E$LightInVoxel[is_sp])
        sp_output[row_nb, ColSMaxLight] <- max(E$LightInVoxel[is_sp])
        sp_output[row_nb, ColSMeanLight] <- mean(E$LightInVoxel[is_sp])
        sp_output[row_nb, ColSMinHeight] <- min(E$Z[is_sp])
        sp_output[row_nb, ColSMaxHeight] <- max(E$Z[is_sp])
        sp_output[row_nb, ColSMeanHeight] <- mean(E$Z[is_sp])
      } else {
        sp_output[row_nb, col_growth_rate] <- NA
        sp_output[row_nb, col_growth_log] <- NA
        sp_output[row_nb, col_birth] <- NA
        sp_output[row_nb, col_death] <- NA
        sp_output[row_nb, col_size] <- NA
        sp_output[row_nb, ColSAverageAge] <- NA
        sp_output[row_nb, ColSMinLight] <- NA
        sp_output[row_nb, ColSMaxLight] <- NA
        sp_output[row_nb, ColSMeanLight] <- NA
        sp_output[row_nb, ColSMinHeight] <- NA
        sp_output[row_nb, ColSMaxHeight] <- NA
        sp_output[row_nb, ColSMeanHeight] <- NA
      }

    } # species loop

    # Community-level output
    MortalityCompetition <- length(which(E$Status == 2))
    MortalityBranchFall <- length(which(E$Status == 3))
    MortalityLight <- length(which(E$Status == 4))
    MortalityNatural <- length(which(E$Status == 5))
    comm_output$timeStep[t] <- InitialTimeStep + t - 1
    comm_output$NumberSpeciesBeginning[t] <- InitialNumberSpecies
    comm_output$NumberSpeciesEnd[t] <- length(unique(E$SpeciesID[E$Status == 1]))
    comm_output$NumberIndividualsBeginning[t] <- nbIndsBeforeDispTotal
    comm_output$NumberIndividualsEnd[t] <- length(which(E$Status == 1))
    comm_output$nb_recruits_matrix[t] <- NumberRecruits
    comm_output$MortalityBranchFall[t] <- MortalityBranchFall
    comm_output$MortalityLight[t] <- MortalityLight
    comm_output$MortalityCompetition[t] <- MortalityCompetition
    comm_output$MortalityNatural[t] <- MortalityNatural
    comm_output$BranchSurfaceIndex[t] <- sum(Microhabitat[, , , 1]) /
      (dimX[1] * dimY[2])
    comm_output$EpiphyteFilling[t] <- sum(E$Mass^(2/3)) /
      SurfaceBiomassScaling / sum(Microhabitat[, , , 1])

    # Command window information
    "--------------------------------------------" |>
      paste_wrap("Species Pool: ", numPool) |>
      paste_wrap("Replicate: ", r) |>
      paste_wrap("Time step: ", InitialTimeStep + t - 1) |>
      paste_wrap("Species Pool: ", numPool) |>
      paste_wrap("Number of individuals: ",
                 comm_output$NumberIndividualsEnd[t]) |>
      paste_wrap("Number of species: ",
                 comm_output$NumberSpeciesEnd[t]) |>
      paste_wrap("Number of recruits: ", NumberRecruits) |>
      paste_wrap("MortalityBranchFall: ", MortalityBranchFall) |>
      paste_wrap("MortalityLight: ", MortalityLight) |>
      paste_wrap("MortalityCompetition: ", MortalityCompetition) |>
      paste_wrap("MortalityNatural: ", MortalityNatural) |>
      paste_wrap("Time: ", format(Sys.time(), "%H:%M:%OS3")) |>
      writeLines()

    # Save Epiphyte matrix for every time step
    ind_output_file <- sub("*.csv$", paste0("_", t, ".csv"), path_to_ind_output)
    write.csv(E[, inds_output_names()], ind_output_file, row.names = FALSE)

    # Save sp_output for every time step
    sp_output_df <- as.data.frame(sp_output)
    names(sp_output_df) <- sp_output_headers
    write.csv(sp_output_df, path_to_sp_output, row.names = FALSE)

    # Save comm_output for every time step (overwrite old one)
    write.csv(comm_output, path_to_comm_output, append = FALSE, row.names = FALSE)

    # Remove dead individuals from Epimatrix
    E <- E[E$Status <= 1, ]

  } # time loop

}
