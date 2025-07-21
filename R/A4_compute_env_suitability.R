options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

setwd("/home/jtrost_ext/MoDVE/R/")

# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("utils.R")

library("doRNG")
library("foreach")
library("doParallel")
# BiocManager::install("rhdf5")
library("rhdf5")


# Parabolic Optimum function
Parabol <- function(a, b, c, x) {
    return((a * x^2) + (b * x) + c)
}


SuitabilityScore <- function (MinEnvVar, MaxEnvVar, OptEnvVar, EnvVar) {

    # Pre-compute denominators
    MaxOptDiff <- MaxEnvVar - OptEnvVar
    OptMinDiff <- OptEnvVar - MinEnvVar

    # Compute suitability only for valid entries
    num <- (MaxEnvVar - EnvVar) / MaxOptDiff
    denom <- (EnvVar - MinEnvVar) / OptMinDiff
    expo  <- OptMinDiff / MaxOptDiff

    suitability <- num * denom^expo

    return(suitability)  # shape: e.g. [50, 50, 60, 100, 2]
}


main <- function() {
    # Detect the number of CPU cores and register the parallel backend
    # Note: detectCores() will detect the total number of cores on a HPC node,
    # so use the $SLURM_NTASKS env var if defined.
    nTasks <- Sys.getenv("SLURM_NTASKS")
    if (nTasks != "") {
        numCores <- strtoi(nTasks)
    }
    else {
        numCores <- detectCores() - 1
    }

    registerDoParallel(numCores)

    # Parse input configuration file
    config <- parse_config()

    DirectoryMicrohabitat <- config$DirectoryMicrohabitat
    DirectorySpeciesPools <- config$DirectorySpeciesPools
    DirectoryOutput <- config$DirectoryOutput
    timeSteps <- config$timeSteps  # Model for timeSteps beginning at the time step given by the initial distribution
    InitialTimeStep <- config$InitialTimeStep
    # Choose species pools to use and number of replicates per species pool
    numSpeciesPools <- config$numSpeciesPools
    # Specifying light response function for growth
    LightResponseFct <- config$LightResponseFct
    Imax <- config$Imax  # Maximum light intensity
    # Get niches to include in the suitability calculation
    LightNicheOpt <- config$LightNicheOpt
    HumNicheOpt <- config$HumNicheOpt
    TempNicheOpt <- config$TempNicheOpt
    WindNicheOpt <- config$WindNicheOpt

    # Load plot dimensions
    dimPlot <- readRDS(file.path(DirectoryMicrohabitat, "dimPlot.rds"))

    # Get microhabitat matrix variables - dynamic handling of selected variables
    MicrohabitatVariableFlags <- config$MicrohabitatVariableFlags
    MhVarNames <- c("TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt", "AverageWeightedAngles",
                    "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
    # Only keep active options
    ActiveOpts <- MhVarNames[as.logical(MicrohabitatVariableFlags)]
    # Assign indices
    Inds <- setNames(seq_along(ActiveOpts), ActiveOpts)
    # Check which light response function to use and get Mh indices for suitability calculation in correct order
    nicheFlags <- c(Light = LightNicheOpt, Hum = HumNicheOpt, Temp = TempNicheOpt, Wind = WindNicheOpt)
    # Filter names where value is 1
    EnvScoreVars <- names(nicheFlags[nicheFlags == 1])

    if (LightResponseFct != "Yan and Hunt") {
        EnvScoreVars <- EnvScoreVars[EnvScoreVars != "Light"]
    }
    MhIdx <- Inds[paste0(EnvScoreVars, "NicheOpt")]
    allEnvVarNames <- c("Light", "Hum", "Temp", "Wind")
    allEnvVarsIdx <- Inds[paste0(allEnvVarNames, "NicheOpt")]

    # Create folder to save the model results
    dir.create(DirectoryOutput, recursive=TRUE)

    # Parallelize across species pools and species
    output <- foreach(numPool=seq(numSpeciesPools),
                       .export=c("ComputeSuitabilityUnscaled", "int_seq", "Parabol", "SuitabilityScore")) %dorng% {

        # -------- Generate individual suitability scores for each species pool --------

        # Load species pool
        SpeciesPoolFileName <- paste0("SpeciesPool", numPool, ".csv")
        if (!file.exists(file.path(DirectorySpeciesPools, SpeciesPoolFileName))) {
            return(NULL)
        }
        SpeciesPool <- read.csv(file.path(DirectorySpeciesPools, SpeciesPoolFileName), sep=",", header=TRUE)
        NSpecies <- nrow(SpeciesPool)

        # Create Save-Directory for each each replicate/initialDistribution
        DirectoryOutputSpeciesPool <- file.path(DirectoryOutput, "EnvSuitability")
        dir.create(DirectoryOutputSpeciesPool, recursive=TRUE)

        print(paste0("Computing suitability scores for species pool ", numPool, "for each variable ..."))

        pb <- txtProgressBar(min = 0, max = (timeSteps + 1), style = 3)

        for (step in seq(0, timeSteps)) {
            t <- InitialTimeStep + step

            savePath <- file.path(DirectoryOutputSpeciesPool,
                                  paste0("ID_SpeciesP_", numPool, "_TimeStep", t, ".h5"))
            if (file.exists(savePath)) {
                h5ds <- h5ls(savePath, all=TRUE)$name
                if ("EnvironmentalSuitabilityScores" %in% h5ds) {
                    # If the file already exists, skip to the next time step
                    setTxtProgressBar(pb, step + 1)
                    next
                } else {
                    # If the file exists but is missing the scores
                    print(paste0("Datasets in ", savePath, ": ", h5ds))
                }
            }

            # Efficiently read once per timestep
            FileNameMicrohabitat <- file.path(DirectoryMicrohabitat, paste0("MicrohabitatMatrix", t, ".rds"))
            if (!file.exists(FileNameMicrohabitat)) {
            stop(paste("Microhabitat file for time step", t, "does not exist:", FileNameMicrohabitat))
            }
            Microhabitat <- readRDS(FileNameMicrohabitat)

            # Scale light
            Microhabitat[, , , Inds["LightNicheOpt"]] <- Imax * Microhabitat[, , , Inds["LightNicheOpt"]]

            SuitabilityScoresT <- array(NA, dim=c(dimPlot, nrow(SpeciesPool), length(allEnvVarsIdx)))

            for (j in seq_along(allEnvVarsIdx)) {
                envVarIdx <- allEnvVarsIdx[j]

                # Extract environmental values in bulk
                envVar <- c(Microhabitat[, , , envVarIdx])
                VarName <- strsplit(names(envVarIdx), split='NicheOpt', fixed=TRUE)[[1]]

                for (i in seq_len(nrow(SpeciesPool))) {

                    # Compute suitability in vectorized form
                    if (VarName == "Light" & LightResponseFct == "Parabolic") {
                        EnvVarSuit <- Parabol(
                            SpeciesPool$LightResponseA[i], SpeciesPool$LightResponseB[i],
                            SpeciesPool$LightResponseC[i], Light
                        )
                    } else {
                        # Use the Yan and Hunt light response function
                        EnvVarSuit <- SuitabilityScore(
                            SpeciesPool[i, paste0("Min", VarName)],
                            SpeciesPool[i, paste0("Max", VarName)],
                            SpeciesPool[i, paste0("Optimum", VarName)],
                            envVar
                        )
                    }
                    EnvVarSuit[is.nan(EnvVarSuit) | is.na(EnvVarSuit)] <- 0  # Set NaN/NA to 0
                    EnvVarSuit[EnvVarSuit < 0] <- 0  # Set negative values to 0
                    SuitabilityScoresT[,,, i, j] <- array(EnvVarSuit, dim=dimPlot)

                    # Print avg. score for this variable
                    avgSuitability <- mean(EnvVarSuit, na.rm=TRUE)
                    print(paste0("Species: ", SpeciesPool$Species[i],
                                 ", Variable: ", VarName,
                                 ", Avg. Suitability: ", round(avgSuitability, 3)))
                }
            }
            # Store suitability scores for the current timestep
            rhdf5::h5createFile(savePath)
            rhdf5::h5write(SuitabilityScoresT, savePath, "EnvironmentalSuitabilityScores")
            setTxtProgressBar(pb, step + 1)
        }
        close(pb)

        # -------- Compute a combined score and scale it --------

        # - 1. Compute the global maximum suitability across all time steps
        globalMaxSuitability <- rep(-Inf, NSpecies)
        activeNiches <- nicheFlags[allEnvVarNames]

        print(paste0("Computing max. suitability scores for species pool ", numPool, "for each species ..."))
        pb <- txtProgressBar(min = 0, max = (timeSteps + 1), style = 3)

        for (step in 0:timeSteps) {
            t <- InitialTimeStep + step
            savePath <- file.path(DirectoryOutputSpeciesPool,
                                  paste0("ID_SpeciesP_", numPool, "_TimeStep", t, ".h5"))
            SuitabilityScoresT <- rhdf5::h5read(savePath, "EnvironmentalSuitabilityScores")

            # Extract the relevant environmental variables
            selectedScores <- SuitabilityScoresT[,,,, activeNiches]

            # Multiply if more than one niche is selected
            if (sum(nicheFlags) > 1) {
                EnvSuitability <- apply(selectedScores, c(1, 2, 3, 4), prod)
            } else {
                EnvSuitability <- selectedScores
            }

            # Get the maximum suitability for this time step for later scaling
            maxThisStep <- apply(EnvSuitability, 4, max, na.rm = TRUE)
            isNewMax <- maxThisStep > globalMaxSuitability
            globalMaxSuitability[isNewMax] <- maxThisStep[isNewMax]

            setTxtProgressBar(pb, step + 1)
        }
        close(pb)
        # guard: if any species were all NA across time, max stays -Inf -> set to NA (or 0)
        globalMaxSuitability[is.infinite(globalMaxSuitability)] <- NA_real_

        # - 2. Recompute suitability scores for each time step and scale them
        print(paste0("Recompute combined scores and scale them for species pool ", numPool, " ..."))
        pb <- txtProgressBar(min = 0, max = (timeSteps + 1), style = 3)

        # Create timestamped directory
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        timestampedDir <- file.path(DirectoryOutputSpeciesPool, timestamp)
        dir.create(timestampedDir, recursive = TRUE)

        for (step in 0:timeSteps) {
            t <- InitialTimeStep + step
            inFile <- file.path(
            DirectoryOutputSpeciesPool,
            paste0("ID_SpeciesP_", numPool, "_TimeStep", t, ".h5")
            )
            outFile <- file.path(
            timestampedDir,
            paste0("ScaledSuitability_", numPool, "_TimeStep", t, ".h5")
            )

            SuitabilityScoresT <- h5read(inFile, "EnvironmentalSuitabilityScores")
            selectedScores <- SuitabilityScoresT[,,,, activeNiches, drop = FALSE]

            if (sum(nicheFlags) > 1) {
                EnvSuitability <- apply(selectedScores, c(1, 2, 3, 4), prod)
            } else {
                EnvSuitability <- selectedScores
            }

            # Scale by species max
            # safe denom: if NA (never observed) -> keep NA; if 0 -> avoid divide-by-zero
            denom <- globalMaxSuitability
            denom[is.na(denom) | denom == 0] <- NA_real_
            scaledSuitability <- sweep(EnvSuitability, 4, denom, "/")

            avgSuitability <- mean(EnvSuitability, na.rm=TRUE)
            print(paste0("Avg. Suitability: ", round(avgSuitability, 3)))
            avgScaledSuitability <- mean(scaledSuitability, na.rm=TRUE)
            print(paste0("Avg. Scaled Suitability: ", round(avgSuitability, 3)))

            # Clamp to [0,1]
            scaledSuitability[is.na(scaledSuitability)] <- 0
            scaledSuitability[is.nan(scaledSuitability)] <- 0
            scaledSuitability[scaledSuitability < 0] <- 0
            scaledSuitability[scaledSuitability > 1] <- 1

            # Save file
            if (file.exists(outFile)) file.remove(outFile)
            h5createFile(outFile)
            h5write(scaledSuitability, outFile, "ScaledSuitabilityScores")

            setTxtProgressBar(pb, step + 1)
        }
        # Store the per-species max used for scaling in the same file
        h5write(globalMaxSuitability, outFile, "GlobalMaxSuitability")
        close(pb)
    }
}

main()