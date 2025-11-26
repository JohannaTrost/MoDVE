options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("utils.R")

library("foreach")
library("doParallel")
library("doRNG")

###############################################################################

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

# Bertalanffy Growth
GrowthRate <- function(MaxMass, Mass, K) {
    return(K * (MaxMass - Mass))
}

# Parabolic Optimum function
Parabol <- function(a, b, c, x) {
    return((a * x^2) + (b * x) + c)
}

dispersal <- function(NumberOfSpecies,
                      E,
                      Microhabitat,
                      SurfaceBiomassScaling,
                      dimPlot,
                      centralPoint,
                      InterceptRecruitment,
                      SlopeRecruitment,
                      prob_disp_matrix,
                      SpeciesPool,
                      max_id) {

    # Store number of individuals at beginning of time step
    nbIndsBeforeDisp <- array(rep(0, NumberOfSpecies))
    for (sp in seq_len(NumberOfSpecies)) {
        nbIndsBeforeDisp[sp] <- length(which(E$SpeciesID == sp & E$Status == 1))
    }
    nbIndsBeforeDispTotal <- length(which(E$Status == 1))
    nbRecruitsPerSpecies <- array(rep(0, NumberOfSpecies))

    # Calculate available surface area per voxel
    avail_sa_matrix <- Microhabitat[, , , 1]
    for (i in seq_len(nrow(E))) {
        sa_needed <- E$Mass[i]^(2/3) / SurfaceBiomassScaling
        avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] <- max(
          0, avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] - sa_needed
          )
    }

    # Initialize potential recruitment dataframe
    unique_species <- unique(E$SpeciesID)
    nbSpecies <- length(unique_species)
    PotentialRecruitment <- data.frame(
      species_index = seq_len(nbSpecies),
      nb_potential_recruits = numeric(nbSpecies)
    )

    # Loop over all species
    for (i in seq_len(nbSpecies)) {

      sp <- unique_species[i]

      # Generate initially empty matrix to store the probabilities for recruitment
      exptd_nb_recruits_matrix <- array(
        rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]),
        dim = c(dimPlot[1], dimPlot[2], dimPlot[3])
      )

      # Matrix containing all mature individuals of one species
      mature_inds <- E[E$SpeciesID == sp & E$Mass >= E$MassAtMaturity, ]

      minLight <- mature_inds$MinLight[1] # same for all species members
      maxLight <- mature_inds$MaxLight[1]

      # Probability matrix for each species:
      # Depending on the position of each mature individual,
      # the total probability for the species is calculated.
      #
      # The second part of the equation accounts for the actual size of the individual
      # in relation to the maximum size for which the recruitment per individual is defined
      for (j in seq_len(nrow(mature_inds))) {

        dist_to_center <- c(
          centralPoint[1] - mature_inds$X[j],
          centralPoint[2] - mature_inds$Y[j],
          centralPoint[3] - mature_inds$Z[j]
          )
        x_coords <- seq(dist_to_center[1] + 1, dist_to_center[1] + dimPlot[1])
        y_coords <- seq(dist_to_center[2] + 1, dist_to_center[2] + dimPlot[2])
        z_coords <- seq(dist_to_center[3] + 1, dist_to_center[3] + dimPlot[3])

        # TODO: none of the following terms are documented in the paper?
        # Mass-dependent fecundity?
        factor1 <- (InterceptRecruitment + SlopeRecruitment * mature_inds$Mass[j]) *
          mature_inds$RecruitmentInvestmentRel[j]
        # ??
        factor2 <- (mature_inds$Mass[j] - mature_inds$MassAtMaturity[j]) /
          (mature_inds$MaximumMass[j] - mature_inds$MassAtMaturity[j])
        # ??
        factor3 <- 1 + (mature_inds$RecruitmentInc[j] * factor2)

        # Dispersal probability * fecundity = expected nb offspring in each xyz
        exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix +
          prob_disp_matrix[x_coords, y_coords, z_coords, sp] * factor1 * factor3
      }

      # Store potential normalized number of recruits
      # We will use this to populate sp_output_mat later
      PotentialRecruitment$nb_potential_recruits[i] <- sum(exptd_nb_recruits_matrix)

      # Matrix containing all voxel for which the light requirements are fulfilled
      pot_hab_matrix <- ifelse(
        Microhabitat[, , , 3] >= minLight & Microhabitat[, , , 3] <= maxLight,
        1, 0
        )

      # Disable unsuitable cells and scale with surface area
      exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix *
        pot_hab_matrix * avail_sa_matrix # TODO: confirm fecundity scales with SA?

      # Calculate number of recruits based on final probability matrix
      nb_recruits_matrix <- array(
        rpois(length(exptd_nb_recruits_matrix), exptd_nb_recruits_matrix),
        dim = dim(exptd_nb_recruits_matrix)
        )

      # Increment recruit counts
      totalNbRecruits <- sum(nb_recruits_matrix)
      nbRecruitsPerSpecies[sp] <- totalNbRecruits

      # Add new recruits to epiphyte matrix
      if (totalNbRecruits > 0) {
        ids <- array_recruits(which(nb_recruits_matrix > 0), dim(nb_recruits_matrix))
        x_recruits <- ids[, 1]
        y_recruits <- ids[, 2]
        z_recruits <- ids[, 3]

        while (totalNbRecruits > length(x_recruits)) {
          # recursively distribute coordinates until counts of remaining
          # unprocessed recruits reaches zero

          tmp_ids <- array_recruits(which(nb_recruits_matrix > 0), dim(nb_recruits_matrix))

          # decrement
          nb_recruits_matrix[tmp_ids] = nb_recruits_matrix[tmp_ids] - 1

          tmp_ids <- array_recruits(which(nb_recruits_matrix > 0), dim(nb_recruits_matrix))

          x_recruits <- append(x_recruits, tmp_ids[, 1])
          y_recruits <- append(y_recruits, tmp_ids[, 2])
          z_recruits <- append(z_recruits, tmp_ids[, 3])
        }

        # Fill Epiphyte matrix
        idx_recruits <- seq(nrow(E) + 1, nrow(E) + length(x_recruits))
        E[idx_recruits, names(SpeciesPool)] <- SpeciesPool[sp, ]
        E$X[idx_recruits] <- x_recruits
        E$Y[idx_recruits] <- y_recruits
        E$Z[idx_recruits] <- z_recruits
        E$Mass[idx_recruits] <- 0  # Initial size
        E$Status[idx_recruits] <- 1  # status 1:alive
        recruits_ids <- seq(max_id + 1, max_id + length(x_recruits))
        E$IndividualID[idx_recruits] <- recruits_ids

        E[is.na(E)] <- 0  # convert all NA to 0 so that the R script matches the Matlab
        # TODO: this is likely to cause bugs

        max_id <- max_id + length(x_recruits)
      } # if any recruits

    } # species loop

    disp_items <- list(
      "nbIndsBeforeDisp" = nbIndsBeforeDisp,
      "nbRecruitsPerSpecies" = nbRecruitsPerSpecies,
      "nbIndsBeforeDispTotal" = nbIndsBeforeDispTotal,
      "E" = E,
      "PotentialRecruitment" = PotentialRecruitment,
      "max_id" = max_id
    )
    return(disp_items)
}


create_pairs <- function(numSpeciesPools, replicatePerSpeciesPool){
    N <- (numSpeciesPools[2] - numSpeciesPools[1] + 1) * replicatePerSpeciesPool
    pairs <- data.frame(matrix(0, nrow=N, ncol=2))
    colnames(pairs) <- c("numPool", "r")

    i <- 1
    for (numPool in int_seq(numSpeciesPools[1], numSpeciesPools[2])) {
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
      numSpeciesPools <- config$numSpeciesPools  # Start and end number of  species pools (if the species pools do not exist, they are automatically skipped)
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
        # If we provide an integer, use it to set the seed
        # otherwise it will be NULL and therefore a random
        # seed will be created.
        seed <- config$seed
        set.seed(seed, kind = "Mersenne-Twister")
    }
    # Save current random state.
    # We need to store this after setting the seed but before calling any
    # functions that generate random numbers.
    writeLines("Storing the random number generator state...")
    save_rng(file.path(DirectoryModelResults, "random_state_seed.RData"))

    # Load plot dimensions
    dimPlot <- readRDS(file.path(DirectoryMicrohabitat, "dimPlot.rds"))

    # Stop criterion: if stop density is exceeded, the simulation ends
    # to cap memory usage
    StopNbInds <- StopDensity * dimPlot[1] * dimPlot[2] / 10000

    # Load TraitRanges (ranges used to create the species pool)
    FileTraitRanges <- file.path(DirectorySpeciesPools, "TraitRanges.csv")
    TraitRanges <- read.table(FileTraitRanges, sep=",", header=FALSE)

    SlopeRecruitment <- TraitRanges[1, 1]
    InterceptRecruitment <- TraitRanges[2, 1]

    # Following the columns in the species matrix refering to a trait or
    # variable. This is handy if the epiphyte matrix changes
    col_sp_id <- 1
    col_nb_inds_begin <- 2
    col_nb_inds_ends <- 3
    col_nb_mature_inds <- 4
    col_nb_recruits <- 5
    ColSPotentialNbRecruits <- 6
    col_nb_dead_branch <- 7
    col_nb_dead_light <- 8
    col_nb_dead_comp <- 9
    col_nb_dead_base <- 10
    ColSNumberPopulationGrowthRate <- 11
    ColSNumberPopulationGrowthRateLog <- 12
    ColSNumberBirthRate <- 13
    ColSNumberDeathRate <- 14
    ColSAverageSize <- 15
    ColSAverageAge <- 16
    ColSMinLight <- 17
    ColSMaxLight <- 18
    ColSMeanLight <- 19
    ColSMinHeight <- 20
    ColSMaxHeight <- 21
    ColSMeanHeight <- 22
    TotalColsSpeciesMatrix <- ColSMeanHeight

    # Headers of matrix
    sp_output_matHeaders <- c("TimeStep", "SpeciesID", "NumberIndividualsBeginning", "NumberIndividualsEnd",
        "NumberMatureIndividuals", "NumberRecruits", "NumberRecruitsPotential", "NumberMortalityBranchFall",
        "NumberMortalityLight", "NumberMortalityCompetition", "NumberMortalityNatural", "PopulationGrowthRate",
        "PopulationGrowthRateLog", "BirthRate", "DeathRate", "AverageMass", "AverageAge", "MinLight", "MaxLight",
        "MeanLight", "MinHeight", "MaxHeight", "MeanHeight")

    if (length(sp_output_matHeaders) != (TotalColsSpeciesMatrix + 1)) {
        stop("Headers of species matrix do not match with number of columns")
    }

    # Headers of matrix
    comm_output_matHeaders <- c("timeStep", "NumberSpeciesBeginning", "NumberSpeciesEnd",
        "NumberIndividualsBeginning", "NumberIndividualsEnd", "nb_recruits_matrix", "MortalityBranchFall",
        "MortalityLight", "MortalityCompetition", "MortalityNatural", "BranchSurfaceIndex", "EpiphyteFilling")

    # Main loop for the community model for each species pool and for each replicate
    # Each loop employs a different random number generator (RNG) stream, resulting in distinct,
    # statistically independent random sequences. These sequences are reproducible across multiple
    # runs, provided the master seed and other input parameters remain unchanged. Importantly, the
    # number of cores does not influence the sequences, only the order in which the loops are
    # processed. Consequently, parallel and serial execution will yield identical results.
    # Internally, the foreach package employs the L'Ecuyer-CMRG RNG algorithm for reliable random
    # number generation, ensuring reproducible results even in parallel computing environments.
    pairs <- create_pairs(numSpeciesPools, replicatePerSpeciesPool)
    shared_objects <- c("calc_prob_disp_matrix", "int_seq", "dispersal", "GrowthRate", "Parabol")
    pair_seq <- seq_len(nrow(pairs))

    output <- foreach::foreach(pair_idx = pair_seq, .export = shared_objects) %dorng% {

        numPool <- pairs$numPool[pair_idx]
        r <- pairs$r[pair_idx]

        # Check if a initial distribution for the species pool exists. If not, move on to the next species pool
        FileNameInitialDistribution <- file.path(DirectoryModelMain, paste("ID_SpeciesP_", numPool, "_Rep_", r, ".csv", sep=""))
        if (!file.exists(FileNameInitialDistribution)) {
            return(NULL)
        }

        # Load species pool
        SpeciesPoolFileName <- paste("SpeciesPool", numPool, ".csv", sep="")
        SpeciesPool <- read.csv(file.path(DirectorySpeciesPools, SpeciesPoolFileName), sep=",", header=TRUE)
        NumberOfSpecies <- nrow(SpeciesPool)  # number of species per 25X25m plot

        # Calculate the probability to disperse in surrounding voxels
        dimX <- dimPlot[1] * 2 + 1
        dimY <- dimPlot[2] * 2 + 1
        dimZ <- dimPlot[3] * 2 + 1
        centralPoint <- c(
          floor(dimX/2) + 1,
          floor(dimY/2) + 1,
          floor(dimZ/2) + 1
          )
        prob_disp_matrix <- calc_prob_disp_matrix(
          centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool
          )

        # Create Save-Directory for each each replicate/initialDistribution
        DirectoryModelResultsRun <- file.path(DirectoryModelResults, paste("ID_SpeciesP_", numPool, "_Rep_", r, sep=""))
        dir.create(DirectoryModelResultsRun, recursive=TRUE)

        # Load initial epiphyte distribution
        # E = epiphytes
        E <- read.csv(FileNameInitialDistribution, sep=",", header=TRUE)
        # Add columns to E for additional info
        E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0
        max_id <- nrow(E)  # to trace individual IDs

        if (MicrohabitatType == 2 || MicrohabitatType == 3) {

            # Microhabitat is loaded once per simulation
            Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, "MicrohabitatMatrix1.rds"))

            # Relative to absolute light values ?mol*m-2*s-1
            Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax

            d1 <- dim(Microhabitat)[1]
            d2 <- dim(Microhabitat)[2]
            d3 <- dim(Microhabitat)[3]
            pot_hab_matrix <- array(
              rep(0, d1 * d2 * d3),
              dim = c(d1, d2, d3)
              )
        }

        # Initialize Matrix where community parameters are saved
        comm_output_mat <- data.frame(matrix(
          0.0,
          nrow = timeSteps,
          ncol = length(comm_output_matHeaders)
        ))
        colnames(comm_output_mat) <- comm_output_matHeaders

        # Initialize matrices where the aggregated info on species level are saved
        sp_output_matSave <- array(
          rep(0, (timeSteps*NumberOfSpecies) * (TotalColsSpeciesMatrix + 1)),
          dim = c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix + 1)
          )

        # Initialize Matrix where species parameters are saved
        sp_output_mat <- array(
          rep(0, (timeSteps*NumberOfSpecies) * TotalColsSpeciesMatrix),
          dim = c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix)
          )

        # Generation loop
        for (t in seq_len(timeSteps)) {

          # Check if the stop criterion is met
          nbIndsAlive <- length(which(E$Status == 1))
          if (nbIndsAlive > StopNbInds) break

          # Update microhabitat if applicable
          if (MicrohabitatType == 1) { # Dynamic forest
            Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, paste("MicrohabitatMatrix", InitialTimeStep + t - 1, ".rds", sep="")))
            Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax
            d1 <- dim(Microhabitat)[1]
            d2 <- dim(Microhabitat)[2]
            d3 <- dim(Microhabitat)[3]
            pot_hab_matrix <- array(rep(0, d1 * d2 * d3), dim=c(d1, d2, d3))
          }

          # Update how many species are alive at beginning of generation
          InitialNumberSpecies <- length(unique(E$SpeciesID[E$Status == 1]))

          # Dispersal
          disp_items <- dispersal(
            NumberOfSpecies,
            E,
            Microhabitat,
            SurfaceBiomassScaling,
            dimPlot,
            centralPoint,
            InterceptRecruitment,
            SlopeRecruitment,
            prob_disp_matrix,
            SpeciesPool,
            max_id
          )

          # Unwrap dispersal output
          nbIndsBeforeDisp <- disp_items$nbIndsBeforeDisp
          nbRecruitsPerSpecies <- disp_items$nbRecruitsPerSpecies
          nbIndsBeforeDispTotal <- disp_items$nbIndsBeforeDispTotal
          PotentialRecruitment <- disp_items$PotentialRecruitment
          E <- disp_items$E
          max_id <- disp_items$max_id

          # Store potential normalized number of recruits in sp_output_mat
          for (ii in seq_len(nrow(PotentialRecruitment))) {
            kk <- PotentialRecruitment$index[ii]
            if (kk != 0) {
              rowID <- (kk-1) * timeSteps + t
              sp_output_mat[rowID, ColSPotentialNbRecruits] <-
                PotentialRecruitment$nb_potential_recruits[ii]
            }
          }
          NumberRecruits <- length(which(E$Status == 1)) - nbIndsBeforeDispTotal

          # TODO:
          # Unclear what this line in the Matlab script is supposed to do.
          # From what I understand, the first column in E ("SpeciesID") takes non-zero values
          # only, so I think that E(:,1)==0 will always be empty.
          # E(E(:,1)==0,:)=[]; %in rare case, some individuals with only zeros are creates, which is wrong. This is to prevent the script to stop.

          # Growth
          for (i in seq_len(nrow(E))) {
            # maybe it is faster if I do not use the if statement => speed testing
            if (E$Status[i] == 1) {
              growth_term <- GrowthRate(E$MaximumMass[i], E$Mass[i], E$GrowthRate[i])
              light_term <- Parabol(E$LightResponseA[i], E$LightResponseB[i], E$LightResponseC[i], Microhabitat[E$X[i], E$Y[i], E$Z[i], 3])
              E$Mass[i] <- E$Mass[i] + max(0, growth_term * light_term)
            }

            # Add info about the voxel to the epiphyte matrix
            E$SurfaceAreaOccupied[i] <- (E$Mass[i]^(2/3)) / SurfaceBiomassScaling
            # TODO: do we really need individual-level copies of these habitat values?
            E$TotalSurfaceInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 1]  # Total surface in voxel
            E$SurfaceLossInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 2]  # Percentage surface loss in this year
            E$LightInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 3]  # Light conditions in voxel
          }

          # Mortality
          for (i in seq_len(nrow(E))) {
            if (E$Status[i] == 1) {

              # The following comparison would fail without the is.nan check,
              # because Microhabitat contains NaNs in some entries and
              # in R a comparison with a NaN returns NA, not a boolean.
              # Note: We call runif repeatedly intentionally. See Issue #16 on Github

              # Branch fall mortality
              if (!is.nan(Microhabitat[E$X[i], E$Y[i], E$Z[i], 2]) &&
                  runif(1, min = 0, max = 1) < Microhabitat[E$X[i], E$Y[i], E$Z[i], 2]
                  ) {
                E$Status[i] <- 3
              } else if (Microhabitat[E$X[i], E$Y[i], E$Z[i], 3] < E$MinLight[i] |
                         Microhabitat[E$X[i], E$Y[i], E$Z[i], 3] > E$MaxLight[i]) {
                # Unsuitable light conditions
                E$Status[i] <- 4
              } else if (MortalityMethod == 0 && runif(1, min=0, max=1) < MortRateRandom) {  # Natural mortality rate
                # Baseline mortality
                E$Status[i] <- 5
              } else if (MortalityMethod == 1 && runif(1, min=0, max=1) < (MortRateMass * (E$Mass[i]^MortRateMassScaling))) {
                # Mass-dependent mortality
                E$Status[i] <- 5
              }
            }
          }

          # Mortality due to competition for space

          # Calculate total surface area occupied by epiphytes per voxel
          TotalSurfaceArePerVoxelOccupied <- array(
            rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]),
            dim = c(dimPlot[1], dimPlot[2], dimPlot[3])
            )
          for (w in seq_len(nrow(E))) {
            if (E$Status[w] == 1) {
              TotalSurfaceArePerVoxelOccupied[E$X[w], E$Y[w], E$Z[w]] <-
                TotalSurfaceArePerVoxelOccupied[E$X[w], E$Y[w], E$Z[w]] +
                E$SurfaceAreaOccupied[w]
              # for each cell, sum sa over all inds in this cell
            }
          }

          # Indices of voxel where total area of epiphytes exceeds the available surface area
          ind_tmp <- array_recruits(
            which(TotalSurfaceArePerVoxelOccupied > Microhabitat[, , , 1]),
            dim(TotalSurfaceArePerVoxelOccupied)
            )
          IndX <- ind_tmp[, 1]
          IndY <- ind_tmp[, 2]
          IndZ <- ind_tmp[, 3]

          for (i in seq_len(length(IndX))) {
            # Get all individuals in thisvoxel
            isInVoxel <- E$X == IndX[i] & E$Y == IndY[i] & E$Z == IndZ[i]
            indsInVoxel <- E[isInVoxel & E$Status == 1, ]

            # Sort them by size (CompetitionMethod=1) or randomly (CompetitionMethod=2)
            if (CompetitionMethod == 1) {
              indsInVoxel <- indsInVoxel[order(indsInVoxel$SurfaceAreaOccupied, decreasing=TRUE), ]
            } else if (CompetitionMethod == 2) {
              indsInVoxel <- indsInVoxel[sample(seq_len(nrow(indsInVoxel))), ]
            }

            totalSurfAreaOcc <- cumsum(indsInVoxel$SurfaceAreaOccupied)
            availSurfArea <- Microhabitat[IndX[i], IndY[i], IndZ[i], 1]
            # The voxel can support this many individuals
            capacity <- length(which(totalSurfAreaOcc <= availSurfArea))

            if (capacity < nrow(indsInVoxel)) {
              # largest n individuals live, the rest die
              seq_beyond_capacity <- int_seq(capacity + 1, nrow(indsInVoxel))
              dead_ids <- indsInVoxel[seq_beyond_capacity, "IndividualID"]
              these_die <- is.element(E$IndividualID, dead_ids)
              E[these_die, "Status"] <- 2
            }
          }

          # Age increment
          E$Age <- E$Age + 1

          # Save number of mortality event
          MortalityCompetition <- length(which(E$Status == 2))
          MortalityBranchFall <- length(which(E$Status == 3))
          MortalityLight <- length(which(E$Status == 4))
          MortalityNatural <- length(which(E$Status == 5))

          # Species-level output
          for (numSpecies in seq_len(NumberOfSpecies)) {

            row_nb <- (numSpecies - 1) * timeSteps + t
            sp_members <- E$SpeciesID == numSpecies
            sp_output_mat[row_nb, col_sp_id] <- numSpecies
            sp_output_mat[row_nb, col_nb_inds_begin] <- nbIndsBeforeDisp[numSpecies]
            sp_output_mat[row_nb, col_nb_inds_ends] <- sum(E$Status == 1 & sp_members, na.rm=TRUE)
            sp_output_mat[row_nb, col_nb_mature_inds] <- sum(E$Status == 1 & sp_members & E$Mass >= E$MassAtMaturity, na.rm=TRUE)
            sp_output_mat[row_nb, col_nb_recruits] <- nbRecruitsPerSpecies[numSpecies]
            sp_output_mat[row_nb, col_nb_dead_branch] <- sum(E$Status == 3 & sp_members, na.rm=TRUE)
            sp_output_mat[row_nb, col_nb_dead_light] <- sum(E$Status == 4 & sp_members, na.rm=TRUE)
            sp_output_mat[row_nb, col_nb_dead_comp] <- sum(E$Status == 2 & sp_members, na.rm=TRUE)
            sp_output_mat[row_nb, col_nb_dead_base] <- sum(E$Status == 5 & sp_members, na.rm=TRUE)

            if (sum(E$Status == 1 & sp_members, na.rm = TRUE) > 0 &&
                nbIndsBeforeDisp[numSpecies] > 0) {
              sp_output_mat[row_nb, ColSNumberPopulationGrowthRate] <- sp_output_mat[row_nb, col_nb_inds_ends] / sp_output_mat[row_nb, col_nb_inds_begin]
              sp_output_mat[row_nb, ColSNumberPopulationGrowthRateLog] <- log(sp_output_mat[row_nb, ColSNumberPopulationGrowthRate])
              sp_output_mat[row_nb, ColSNumberBirthRate] <- nbRecruitsPerSpecies[numSpecies] / nbIndsBeforeDisp[numSpecies]
              sp_output_mat[row_nb, ColSNumberDeathRate] <- (sum(E$Status == 3 & sp_members, na.rm=TRUE) + sum(E$Status == 4 & sp_members, na.rm=TRUE) + sum(E$Status == 2 & sp_members, na.rm=TRUE) + sum(E$Status == 5 & sp_members, na.rm=TRUE)) / nbIndsBeforeDisp[numSpecies]
              sp_output_mat[row_nb, ColSAverageSize] <- mean(E$Mass[sp_members])
              sp_output_mat[row_nb, ColSAverageAge] <- mean(E$Age[sp_members])
              sp_output_mat[row_nb, ColSMinLight] <- min(E$LightInVoxel[sp_members])
              sp_output_mat[row_nb, ColSMaxLight] <- max(E$LightInVoxel[sp_members])
              sp_output_mat[row_nb, ColSMeanLight] <- mean(E$LightInVoxel[sp_members])
              sp_output_mat[row_nb, ColSMinHeight] <- min(E$Z[sp_members])
              sp_output_mat[row_nb, ColSMaxHeight] <- max(E$Z[sp_members])
              sp_output_mat[row_nb, ColSMeanHeight] <- mean(E$Z[sp_members])
            } else {
              sp_output_mat[row_nb, ColSNumberPopulationGrowthRate] <- NaN
              sp_output_mat[row_nb, ColSNumberPopulationGrowthRateLog] <- NaN
              sp_output_mat[row_nb, ColSNumberBirthRate] <- NaN
              sp_output_mat[row_nb, ColSNumberDeathRate] <- NaN
              sp_output_mat[row_nb, ColSAverageSize] <- NaN
              sp_output_mat[row_nb, ColSAverageAge] <- NaN
              sp_output_mat[row_nb, ColSMinLight] <- NaN
              sp_output_mat[row_nb, ColSMaxLight] <- NaN
              sp_output_mat[row_nb, ColSMeanLight] <- NaN
              sp_output_mat[row_nb, ColSMinHeight] <- NaN
              sp_output_mat[row_nb, ColSMaxHeight] <- NaN
              sp_output_mat[row_nb, ColSMeanHeight] <- NaN
            }

            sp_output_matSave[row_nb, 1] <- InitialTimeStep + t - 1
            sp_output_matSave[row_nb, int_seq(2, TotalColsSpeciesMatrix + 1)] <- sp_output_mat[row_nb, ]

          } # species loop

          # Community-level output
          comm_output_mat$timeStep[t] <- InitialTimeStep + t - 1
          comm_output_mat$NumberSpeciesBeginning[t] <- InitialNumberSpecies
          comm_output_mat$NumberSpeciesEnd[t] <- length(unique(E$SpeciesID[E$Status == 1]))
          comm_output_mat$NumberIndividualsBeginning[t] <- nbIndsBeforeDispTotal
          comm_output_mat$NumberIndividualsEnd[t] <- length(which(E$Status == 1))
          comm_output_mat$nb_recruits_matrix[t] <- NumberRecruits
          comm_output_mat$MortalityBranchFall[t] <- MortalityBranchFall
          comm_output_mat$MortalityLight[t] <- MortalityLight
          comm_output_mat$MortalityCompetition[t] <- MortalityCompetition
          comm_output_mat$MortalityNatural[t] <- MortalityNatural
          comm_output_mat$BranchSurfaceIndex[t] <- sum(Microhabitat[, , , 1]) /
            (dimPlot[1] * dimPlot[2])
          comm_output_mat$EpiphyteFilling[t] <- sum(E$Mass^(2/3)) /
            SurfaceBiomassScaling / sum(Microhabitat[, , , 1])

          # Command window information
          info <- "--------------------------------------------"
          info <- paste(info, paste("Species Pool: ", numPool, sep=""), sep="\n")
          info <- paste(info, paste("Replicate: ", r, sep=""), sep="\n")
          info <- paste(info, paste("Time step: ", InitialTimeStep + t - 1, sep=""), sep="\n")
          info <- paste(info, paste("Number of individuals: ", comm_output_mat$NumberIndividualsEnd[t], sep=""), sep="\n")
          info <- paste(info, paste("Number of species: ", comm_output_mat$NumberSpeciesEnd[t], sep=""), sep="\n")
          info <- paste(info, paste("Number of recruits: ", NumberRecruits, sep=""), sep="\n")
          info <- paste(info, paste("MortalityBranchFall: ", MortalityBranchFall, sep=""), sep="\n")
          info <- paste(info, paste("MortalityLight: ", MortalityLight, sep=""), sep="\n")
          info <- paste(info, paste("MortalityCompetition: ", MortalityCompetition, sep=""), sep="\n")
          info <- paste(info, paste("MortalityNatural: ", MortalityNatural, sep=""), sep="\n")
          info <- paste(info, paste("Time: ", format(Sys.time(), "%H:%M:%OS3"), sep=""), sep="\n")
          writeLines(info)

          # Saving
          # Save Epiphyte matrix for every time step
          ColumsToSave <- c("SpeciesID", "IndividualID", "Status", "Mass", "Age", "X", "Y", "Z", "TotalSurfaceInVoxel", "SurfaceLossInVoxel", "LightInVoxel")
          write.csv(E[, ColumsToSave], file.path(DirectoryModelResultsRun, paste("IndividualMatrixTimeStep", InitialTimeStep + t - 1, ".csv", sep="")), row.names=FALSE)

          # Save sp_output_mat for every time step
          sp_output_df <- as.data.frame(sp_output_matSave)
          names(sp_output_df) <- sp_output_matHeaders
          write.csv(sp_output_df, file.path(DirectoryModelResultsRun, "SpeciesSummary.csv"), row.names=FALSE)

          # Save comm_output_mat for every time step (overwrite old one)
          write.csv(comm_output_mat, file.path(DirectoryModelResultsRun, "CommunitySummary.csv"), row.names=FALSE)

          ###############################################################################

          # Remove dead individuals from Epimatrix
          E <- E[E$Status <= 1, ]  # Remove rows where Status > 1

        } # time loop

        return(NULL)

    } # foreach speciesPool x replicate
}


main()
