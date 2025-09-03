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
    nTasks <- Sys.getenv("SLURM_CPUS_PER_TASK")
    if (nTasks != "") {
        numCores <- strtoi(nTasks)
    }
    else {
        numCores <- detectCores() - 1
    }

    registerDoParallel(numCores)

    writeLines(paste("Using", numCores, "cores for parallel processing."))

    # Parse input configuration file
    args <- commandArgs(trailingOnly = TRUE)
    configFile <- args[1]
    singleStep <- as.integer(args[2])  # the timestep you want to compute

    config <- read.config(configFile)

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
    #allEnvVarNames <- c("Light", "Hum", "Temp", "Wind")
    allEnvVarsIdx <- Inds[paste0(EnvScoreVars, "NicheOpt")]

    # Create folder to save the model results
    dir.create(DirectoryOutput, recursive=TRUE)

    # Create Save-Directory
    DirectoryOutputSpeciesPool <- file.path(DirectoryOutput, "EnvSuitability")
    dir.create(DirectoryOutputSpeciesPool, recursive=TRUE)

    # Create timestamped directory
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    timestampedDir <- file.path(DirectoryOutputSpeciesPool, timestamp)
    dir.create(timestampedDir, recursive = TRUE)

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

        if (!is.na(singleStep)) {

            print(paste0("Computing suitability scores for species pool ", numPool, "for each variable ..."))

            t <- singleStep
            print(paste0("Time step", t))

            savePath <- file.path(DirectoryOutputSpeciesPool,
                                  paste0("ID_SpeciesP_", numPool, "_TimeStep", t, ".h5"))

            # Efficiently read microhabitat for this timestep
            FileNameMicrohabitat <- file.path(DirectoryMicrohabitat, paste0("MicrohabitatMatrix", t, ".rds"))
            if (!file.exists(FileNameMicrohabitat)) {
                stop(paste("Microhabitat file for time step", t, "does not exist:", FileNameMicrohabitat))
            }
            Microhabitat <- readRDS(FileNameMicrohabitat)

            # Scale light
            Microhabitat[, , , Inds["LightNicheOpt"]] <- Imax * Microhabitat[, , , Inds["LightNicheOpt"]]

            SuitabilityScoresT <- array(NA, dim=c(dimPlot, nrow(SpeciesPool), length(allEnvVarsIdx)))

            for (j in seq_along(EnvScoreVars)) {
                envVarIdx <- allEnvVarsIdx[j]
                envVar <- c(Microhabitat[, , , envVarIdx])
                VarName <- strsplit(names(envVarIdx), split='NicheOpt', fixed=TRUE)[[1]]

                for (i in seq_len(nrow(SpeciesPool))) {
                    if (VarName == "Light" & LightResponseFct == "Parabolic") {
                        EnvVarSuit <- Parabol(
                            SpeciesPool$LightResponseA[i], SpeciesPool$LightResponseB[i],
                            SpeciesPool$LightResponseC[i], envVar
                        )
                    } else {
                        EnvVarSuit <- SuitabilityScore(
                            SpeciesPool[i, paste0("Min", VarName)],
                            SpeciesPool[i, paste0("Max", VarName)],
                            SpeciesPool[i, paste0("Optimum", VarName)],
                            envVar
                        )
                    }
                    EnvVarSuit[is.nan(EnvVarSuit) | is.na(EnvVarSuit)] <- 0
                    EnvVarSuit[EnvVarSuit < 0] <- 0
                    SuitabilityScoresT[,,, i, j] <- array(EnvVarSuit, dim=dimPlot)

                    avgSuitability <- mean(EnvVarSuit, na.rm=TRUE)
                    print(paste0("Species: ", SpeciesPool$Species[i],
                                 ", Variable: ", VarName,
                                 ", Avg. Suitability: ", round(avgSuitability, 3)))
                }
            }

            # Save unscaled suitability (no scaling done in single-step mode)
            if (file.exists(savePath)) file.remove(savePath)
            rhdf5::h5createFile(savePath)
            rhdf5::h5write(SuitabilityScoresT, savePath, "EnvironmentalSuitabilityScores")

        } else {

            # -------- Compute a combined score and scale it --------

            # - 1. Compute the global maximum suitability across all time steps
            globalMaxSuitability <- rep(-Inf, NSpecies)
            activeNiches <- nicheFlags[EnvScoreVars]

            print(paste0("Computing max. suitability scores for species pool ", numPool, " for each species ..."))

            # Parallel loop across time steps
            maxList <- foreach(step = 0:timeSteps, .combine = 'cbind', .packages = "rhdf5") %dopar% {
                t <- InitialTimeStep + step

                writeLines(paste("Processing time step", t))

                savePath <- file.path(DirectoryOutputSpeciesPool,
                                      paste0("ID_SpeciesP_", numPool, "_TimeStep", t, ".h5"))
                SuitabilityScoresT <- rhdf5::h5read(savePath, "EnvironmentalSuitabilityScores")

                # Extract the relevant environmental variables
                selectedScores <- SuitabilityScoresT[,,,, activeNiches, drop = FALSE]

                # Multiply if more than one niche is selected
                if (sum(nicheFlags) > 1) {
                    EnvSuitability <- apply(selectedScores, c(1, 2, 3, 4), prod)
                } else {
                    EnvSuitability <- selectedScores
                }

                # Get max per species
                maxThisStep <- apply(EnvSuitability, 4, max, na.rm = TRUE)
                return(maxThisStep)
            }

            # Collapse across steps → take global max
            globalMaxSuitability <- apply(maxList, 1, max, na.rm = TRUE)

            # Guard: if any species were all NA across time, set to NA
            globalMaxSuitability[is.infinite(globalMaxSuitability)] <- NA_real_

            # Store the per-species max used for scaling in the last file
            maxSuitFile <- file.path(
                timestampedDir,
                paste0("GlobalMaxSuitability_", numPool, ".h5")
            )
            tryCatch({
                rhdf5::h5createFile(maxSuitFile)
                rhdf5::h5write(globalMaxSuitability, maxSuitFile, "GlobalMaxSuitability")
            }, error = function(e) {
                message("❌ Failed to save: ", maxSuitFile)
                message("   Error: ", conditionMessage(e))
                return(NA)  # mark this iteration as failed
            })

            # - 2. Recompute suitability scores for each time step and scale them
            print(paste0("Recompute combined scores and scale them for species pool ", numPool, " ..."))

            dummy <- foreach(step = 0:timeSteps, .packages = "rhdf5") %dopar% {
                t <- InitialTimeStep + step
                inFile <- file.path(
                    DirectoryOutputSpeciesPool,
                    paste0("ID_SpeciesP_", numPool, "_TimeStep", t, ".h5")
                )
                outFile <- file.path(
                    timestampedDir,
                    paste0("ScaledSuitability_", numPool, "_TimeStep", t, ".h5")
                )

                SuitabilityScoresT <- rhdf5::h5read(inFile, "EnvironmentalSuitabilityScores")
                selectedScores <- SuitabilityScoresT[,,,, activeNiches, drop = FALSE]

                if (sum(nicheFlags) > 1) {
                    EnvSuitability <- apply(selectedScores, c(1, 2, 3, 4), prod)
                } else {
                    EnvSuitability <- selectedScores
                }

                # Scale by species max
                denom <- globalMaxSuitability
                denom[is.na(denom) | denom == 0] <- NA_real_
                scaledSuitability <- sweep(EnvSuitability, 4, denom, "/")

                avgSuitability <- mean(EnvSuitability, na.rm=TRUE)
                cat("Step", t, ": Avg. Suitability =", round(avgSuitability, 3), "\n")

                avgScaledSuitability <- mean(scaledSuitability, na.rm=TRUE)
                cat("Step", t, ": Avg. Scaled Suitability =", round(avgScaledSuitability, 3), "\n")

                # Clamp to [0,1]
                scaledSuitability[is.na(scaledSuitability)] <- 0
                scaledSuitability[is.nan(scaledSuitability)] <- 0
                scaledSuitability[scaledSuitability < 0] <- 0
                scaledSuitability[scaledSuitability > 1] <- 1

                # --- Safe file write ---
                tryCatch({
                    if (file.exists(outFile)) file.remove(outFile)
                    rhdf5::h5createFile(outFile)
                    rhdf5::h5write(as.array(scaledSuitability), outFile, "ScaledSuitabilityScores")
                }, error = function(e) {
                    message("❌ Failed to save: ", outFile)
                    message("   Error: ", conditionMessage(e))
                    return(NA)  # mark this iteration as failed
                })

                # return file name for debugging if needed
                return(outFile)
            }
        }
    }
}

main()
