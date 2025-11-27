options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("R/utils.R")

find_central_point <- function(dims) {
  dimX <- dims[1] * 2 + 1
  dimY <- dims[2] * 2 + 1
  dimZ <- dims[3] * 2 + 1
  return(c(
    floor(dimX/2) + 1,
    floor(dimY/2) + 1,
    floor(dimZ/2) + 1
  ))
}

calc_prob_disp_matrix <- function(centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool) {

  # Calculate distance to central point
    DistanceMatrix <- array(
      rep(0, dimX * dimY * dimZ),
      dim = c(dimX, dimY, dimZ)
      )
    for (i in seq_len(dimX)) {
        for (j in seq_len(dimY)) {
            for (k in seq_len(dimZ)) {
                x1 <- c(i, j, k)
                x2 <- c(centralPoint[1], centralPoint[2], centralPoint[3])
                DistanceMatrix[i, j, k] <- sqrt(sum((x1 - x2)^2))  # call to pdist() in the matlab script
            }
        }
    }

    # Get probabilities to disperse in each voxel
    ProbabilityMatrix <- prob_disp_matrix <- array(
      rep(0, dimX * dimY * dimZ * NumberOfSpecies),
      dim = c(dimX, dimY, dimZ, NumberOfSpecies)
      )

    for (i in seq_len(NumberOfSpecies)) {

        # Negative exponential
        exponentE <- SpeciesPool$DispersalKernel[i]
        ProbabilityMatrix[, , , i] <- exp(-DistanceMatrix * exponentE)

        # Dispersal asymmetry (probability to disperse downwards > upwards)
        dispersalAsymmetry <- SpeciesPool$DispersalKernelAsymmetry[i]
        z_seq_up <- int_seq(from = centralPoint[3], to = dimZ, by = 1)
        z_seq_down <- int_seq(from = 1, to = centralPoint[3] - 1, by = 1)
        ProbabilityMatrix[, , z_seq_up, i] <- ProbabilityMatrix[, , z_seq_up, i] *
          2 * (1 - dispersalAsymmetry)
        ProbabilityMatrix[, , z_seq_down, i] <- ProbabilityMatrix[, , z_seq_down, i] *
          2 * dispersalAsymmetry

        # Normalize
        prob_disp_matrix[, , , i] <- ProbabilityMatrix[, , , i] /
          sum(ProbabilityMatrix[, , , i])
    }
    return(prob_disp_matrix)
}

# Functions used in the model

create_pairs <- function(nb_spPools, replicatePerSpeciesPool){
    N <- (nb_spPools[2] - nb_spPools[1] + 1) * replicatePerSpeciesPool
    pairs <- data.frame(matrix(0, nrow=N, ncol=2))
    colnames(pairs) <- c("numPool", "r")

    i <- 1
    for (numPool in int_seq(nb_spPools[1], nb_spPools[2])) {
        for (r in seq_len(replicatePerSpeciesPool)) {
            pairs$numPool[i] <- numPool
            pairs$r[i] <- r
            i <- i + 1
        }
    }

    return(pairs)
}

save_rng <- function(savefile) {
    if (exists(".Random.seed")) {
        oldseed <- get(".Random.seed", .GlobalEnv)
    } else {
        stop("You need to call set.seed() first.")
    }
    oldRNGkind <- RNGkind()
    save("oldseed", "oldRNGkind", file=savefile)
}

restore_rng <- function(savefile) {
    load(savefile)
    do.call("RNGkind", as.list(oldRNGkind))
    assign(".Random.seed", oldseed, .GlobalEnv)
}

main <- function() {

    # Detect the number of CPU cores and register the parallel backend
    # Note: detectCores() will detect the total number of cores on a HPC node,
    # so use the $SLURM_NTASKS env var if defined.
    nTasks <- Sys.getenv("SLURM_NTASKS")
    if (nTasks != "") numCores <- strtoi(nTasks)
    else numCores <- parallel::detectCores()
    doParallel::registerDoParallel(numCores)

    # Parse input configuration file
    config <- parse_config("tests/config_a4.toml")

    # Parameters that need to be specified/checked before running this script
    {
      # Input directories
      DirectoryMicrohabitat <- config$DirectoryMicrohabitat
      DirectorySpeciesPools <- config$DirectorySpeciesPools
      DirectoryModelMain <- config$DirectoryModelMain

      # Output directory
      DirectoryModelResults <- config$DirectoryModelResults

      # Define which type of forest the microhabitat belongs to.
      # 1: dynamic forest, 2: static forest, 3: uniform forest
      MicrohabitatType <- config$MicrohabitatType

      # Model parameters

      # Model for timeSteps beginning at the time step given by the initial distribution
      timeSteps <- config$timeSteps

      # Density of individuals per ha at which to stop the simulation of the community and
      # move to the next replicate (to prevent exploding communities)
      StopDensity <- config$StopCriterionHa  # Individuals per ha

      # Choose species pools to use and number of replicates per species pool
      nb_spPools <- config$nb_spPools  # Start and end number of  species pools (if the species pools do not exist, they are automatically skipped)
      replicatePerSpeciesPool <- config$replicatePerSpeciesPool  # Number of replicates per species pool  (if the replicates do not exist, they are automatically skipped)

      SurfaceBiomassScaling <- config$SurfaceBiomassScaling  # cm^2 per m^2
      Imax <- config$Imax  # maximum light above canopy

      # Competition Methods; defines which individuals are removed in voxels which
      # are entirely filled.
      # 1:size (small individuals are outcompetet by larger ones);
      # 2:random competition
      CompetitionMethod <- config$CompetitionMethod

      # Mortality method (complete random or scaling with mass according to metabolic theory);
      # 0: random mortality; 1: scaling with mass to the exponent -1/4
      MortalityMethod <- config$MortalityMethod
      MortRateRandom <- config$MortRateRandom
      MortRateMass <- config$MortRateMass
      MortRateMassScaling <- config$MortRateMassScaling  # widely used scaling factor

      InitialTimeStep <- config$InitialTimeStep  # Time step for which the Initial distribution is generated in A3
    }

    # Create folder to save the model results
    dir.create(DirectoryModelResults, recursive=TRUE)

    # RNG seed
    random_state_file <- config$RandomState
    if (!is.null(random_state_file) && file.exists(random_state_file)) {
      writeLines("Loading previous random number generator state from ", random_state_file)
      restore_rng(random_state_file)
    } else {
      # If we provide an integer, use it to set the seed otherwise it will
      # be NULL and therefore a random seed will be created.
      seed <- config$seed
      set.seed(seed, kind = "Mersenne-Twister")
    }
    # Save current random state.
    # We need to store this after setting the seed but before calling any
    # functions that generate random numbers.
    writeLines("Storing the random number generator state...")
    save_rng(file.path(DirectoryModelResults, "random_state_seed.RData"))

    # Load TraitRanges (ranges used to create the species pool)
    FileTraitRanges <- file.path(DirectorySpeciesPools, "TraitRanges.csv")
    TraitRanges <- read.table(FileTraitRanges, sep=",", header=FALSE)

    # Global parameters of the simulation
    SlopeRecruitment <- TraitRanges[1, 1]
    InterceptRecruitment <- TraitRanges[2, 1]
    # TODO: add these parameters to species pool so we don't load another file just for them

    # Following the columns in the species matrix refering to a trait or
    # variable. This is handy if the epiphyte matrix changes
    {
      col_sp_id <- 1
      col_nb_inds_begin <- 2
      col_nb_inds_end <- 3
      col_nb_mature_inds <- 4
      col_nb_rec <- 5
      col_nb_rec <- 6
      col_nb_dead_branch <- 7
      col_nb_dead_light <- 8
      col_nb_dead_comp <- 9
      col_nb_dead_base <- 10
      col_growth_rate <- 11
      col_growth_log <- 12
      col_birth <- 13
      col_death <- 14
      col_size <- 15
      ColSAverageAge <- 16
      ColSMinLight <- 17
      ColSMaxLight <- 18
      ColSMeanLight <- 19
      ColSMinHeight <- 20
      ColSMaxHeight <- 21
      ColSMeanHeight <- 22
      nb_cols_sp_output <- ColSMeanHeight
    }
    # Headers of matrix
    sp_outputHeaders <- species_output_names()
    if (length(sp_outputHeaders) != (nb_cols_sp_output + 1)) {
        stop("Headers of species matrix do not match with number of columns")
    }

    # Headers of matrix
    comm_outputHeaders <- comm_output_names()

    # Main loop for the community model for each species pool and for each replicate
    # Each loop employs a different random number generator (RNG) stream, resulting in distinct,
    # statistically independent random sequences. These sequences are reproducible across multiple
    # runs, provided the master seed and other input parameters remain unchanged. Importantly, the
    # number of cores does not influence the sequences, only the order in which the loops are
    # processed. Consequently, parallel and serial execution will yield identical results.
    # Internally, the foreach package employs the L'Ecuyer-CMRG RNG algorithm for reliable random
    # number generation, ensuring reproducible results even in parallel computing environments.
    pairs <- create_pairs(nb_spPools, replicatePerSpeciesPool)
    shared_objects <- c("calc_prob_disp_matrix", "int_seq", "resolveDispersal")
    pair_seq <- seq_len(nrow(pairs))

    output <- foreach::foreach(pair_idx = pair_seq, .export = shared_objects) %dorng% {

        numPool <- pairs$numPool[pair_idx]
        r <- pairs$r[pair_idx]

        # Check if a initial distribution for the species pool exists. If not, move on to the next species pool
        FileNameInitialDistribution <- file.path(DirectoryModelMain, paste("ID_SpeciesP_", numPool, "_Rep_", r, ".csv", sep=""))
        if (!file.exists(FileNameInitialDistribution)) return(NULL)

        # Load species pool
        SpeciesPoolFileName <- paste("SpeciesPool", numPool, ".csv", sep="")
        SpeciesPool <- read.csv(file.path(DirectorySpeciesPools, SpeciesPoolFileName), sep=",", header=TRUE)
        NumberOfSpecies <- nrow(SpeciesPool)  # number of species per 25X25m plot

        # Load microhabitat
        # if static then microhabitat is loaded once per simulation
        # else (MicroHabitatType == 1) first habitat used to read dimensions
        Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, "MicrohabitatMatrix1.rds"))

        # Relative to absolute light values ?mol*m-2*s-1
        Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax
        dimX <- dim(Microhabitat)[1]
        dimY <- dim(Microhabitat)[2]
        dimZ <- dim(Microhabitat)[3]

        # Stop criterion: if stop density is exceeded, the simulation ends
        # to cap memory usage
        StopNbInds <- StopDensity * dimX * dimY / 10000

        # Calculate the probability to disperse in surrounding voxels
        centralPoint <- find_central_point(dims)
        prob_disp_matrix <- calc_prob_disp_matrix(
          centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool
        )

        # Create Save-Directory for each each replicate/initialDistribution
        DirectoryModelResultsRun <- file.path(DirectoryModelResults, paste("ID_SpeciesP_", numPool, "_Rep_", r, sep=""))
        dir.create(DirectoryModelResultsRun, recursive=TRUE)

        # Load initial epiphyte distribution
        # E = epiphytes
        E <- read.csv(FileNameInitialDistribution, sep = ",", header = TRUE)
        # Add columns to E for additional info
        E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0
        max_id <- nrow(E)  # to trace individual IDs

        # Initialize Matrix where community parameters are saved
        comm_output <- data.frame(matrix(
          0.0, nrow = timeSteps, ncol = length(comm_outputHeaders)
        ))
        colnames(comm_output) <- comm_outputHeaders

        # Initialize matrices where the aggregated info on species level are saved
        sp_outputSave <- array(
          rep(0, (timeSteps*NumberOfSpecies) * (nb_cols_sp_output + 1)),
          dim = c(timeSteps*NumberOfSpecies, nb_cols_sp_output + 1)
          )

        # Initialize Matrix where species parameters are saved
        sp_output <- array(
          rep(0, (timeSteps*NumberOfSpecies) * nb_cols_sp_output),
          dim = c(timeSteps*NumberOfSpecies, nb_cols_sp_output)
          )

        # Generation loop
        for (t in seq_len(timeSteps)) {

          # Check if the stop criterion is met
          nbIndsAlive <- length(which(E$Status == 1))
          if (nbIndsAlive > StopNbInds) break

          # Update microhabitat if applicable
          if (MicrohabitatType == 1 && t > 1) { # Dynamic forest
            Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, paste("MicrohabitatMatrix", InitialTimeStep + t - 1, ".rds", sep="")))
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

            sp_outputSave[row_nb, 1] <- InitialTimeStep + t - 1
            sp_outputSave[row_nb, int_seq(2, nb_cols_sp_output + 1)] <- sp_output[row_nb, ]

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
          write.csv(
            E[, inds_output_names()],
            file.path(DirectoryModelResultsRun, paste("IndividualMatrixTimeStep", InitialTimeStep + t - 1, ".csv", sep="")),
            row.names = FALSE
            )

          # Save sp_output for every time step
          sp_output_df <- as.data.frame(sp_outputSave)
          names(sp_output_df) <- sp_outputHeaders
          write.csv(
            sp_output_df,
            file.path(DirectoryModelResultsRun, "SpeciesSummary.csv"),
            row.names = FALSE
            )

          # Save comm_output for every time step (overwrite old one)
          write.csv(
            comm_output,
            file.path(DirectoryModelResultsRun, "CommunitySummary.csv"),
            append = FALSE,
            row.names =FALSE
            )

          # Remove dead individuals from Epimatrix
          E <- E[E$Status <= 1, ]

        } # time loop

        return(NULL)

    } # foreach speciesPool x replicate
}
