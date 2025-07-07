options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

#setwd("/home/jtrost_ext/MoDVE/R/")

# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("utils.R")

library("doRNG")
library("foreach")
library("doParallel")
# BiocManager::install("rhdf5")
library("rhdf5")

ComputeSuitabilityUnscaled <- function (timeSteps,
                                   DirectoryMicrohabitat,
                                   SpeciesPool,
                                   InitialTimeStep,
                                   Imax,
                                   LightResponseFct,
                                   Inds,
                                   EnvScoreVars,
                                   DirectoryOutput,
                                   spacialDim) {
    NSpecies <- nrow(SpeciesPool)
    globalMaxSuitability <- rep(-Inf, NSpecies)

    LightIdx <- Inds["LightNicheOpt"]

    # Extract environmental variables from the species pool
    MinEnvVar <- as.matrix(SpeciesPool[paste0("Min", EnvScoreVars)])
    MaxEnvVar <- as.matrix(SpeciesPool[paste0("Max", EnvScoreVars)])
    OptEnvVar <- as.matrix(SpeciesPool[paste0("Optimum", EnvScoreVars)])

    # Output file
    envSuitMatName <- paste0("UnscaledEnvSuitability_t", InitialTimeStep, "-t",
                             InitialTimeStep + timeSteps - 1)
    scaledSuitMatName <- paste0("ScaledEnvSuitability_t", InitialTimeStep, "-t",
                             InitialTimeStep + timeSteps - 1)
    envSuitPath <- file.path(DirectoryOutput, paste0(envSuitMatName, ".h5"))
    scaledSuitPath <- file.path(DirectoryOutput, paste0(scaledSuitMatName, ".h5"))
    maxSuitPath <- file.path(DirectoryOutput,
                             paste0("MaxSuitability_t", InitialTimeStep, "-t",
                                    InitialTimeStep + timeSteps - 1, ".txt"))

    rhdf5::h5createFile(envSuitPath)
    rhdf5::h5createFile(scaledSuitPath)

    fullEnvSuitability <- array(0, dim = c(spacialDim, NSpecies, timeSteps))

    for (t in seq_len(timeSteps)) {

        start <- Sys.time()

        MHPath <- file.path(DirectoryMicrohabitat,
                            paste0("MicrohabitatMatrix", InitialTimeStep + t - 1, ".rds"))
        Microhabitat <- readRDS(MHPath)
        Microhabitat[, , , LightIdx] <- Microhabitat[, , , LightIdx] * Imax

        # -- Compute the score for each environmental variable
        spatialDim <- dim(Microhabitat)
        # Initialize the array for environmental suitability scores with zeros
        EnvSuitabilityVars <- array(
          rep(0, spatialDim[1] * spatialDim[2] * spatialDim[3] * NSpecies * 4),
          dim=c(spatialDim[1], spatialDim[2], spatialDim[3], NSpecies, 4))
        # Compute score for either Hum, Temp, Wind, and Light or only Hum, Temp, and Wind
        EnvSuitabilityVars[, , , , seq_along(EnvScoreVars)] <- SuitabilityScore(
          MinEnvVar, MaxEnvVar, OptEnvVar,
          Microhabitat[ , , , Inds[paste0(EnvScoreVars, "NicheOpt")]])

        if (LightResponseFct == "Parabolic") {
            EnvSuitabilityVars[, , , , 4] <- Parabol(
              SpeciesPool$LightResponseA, SpeciesPool$LightResponseB,
              SpeciesPool$LightResponseC, Microhabitat[ , , , 3])
        }

        # Combine the suitability probabilities for all environmental variables
        EnvSuitability <- apply(EnvSuitabilityVars, c(1, 2, 3, 4), prod)

        # Populate the full environmental suitability array
        fullEnvSuitability[,,,, t] <- EnvSuitability

        # Get the maximum suitability for this time step for later scaling
        maxThisStep <- apply(EnvSuitability, c(4), max, na.rm = TRUE)
        isNewMax <- maxThisStep > globalMaxSuitability
        globalMaxSuitability[isNewMax] <- maxThisStep[isNewMax]

        # Save to disk
        if (t %% 10 == 0 || t == timeSteps) {
            # Save the env suitability for this time step
            rhdf5::h5write(fullEnvSuitability, envSuitPath, envSuitMatName)
            print(paste0("Saving env. suitability up to time step: ", InitialTimeStep + t - 1))
        }

        end <- Sys.time()
        print(paste("Time step", InitialTimeStep + t - 1, "completed in",
                    round(difftime(end, start, units = "secs"), 2), "seconds."))
    }

    # Save the global maximum suitability for scaling to a txt file
    write.table(globalMaxSuitability, maxSuitPath, row.names = FALSE, col.names = FALSE)

    # Scale the suitability scores by the global maximum suitability
    for (s in seq_len(NSpecies)) {
        # Scale the suitability scores for each species
        fullEnvSuitability[,,,s,] <- fullEnvSuitability[,,,s,] / globalMaxSuitability[s]
    }
    # Save the env suitability for this time step
    rhdf5::h5write(fullEnvSuitability, scaledSuitPath, scaledSuitMatName)
    print(paste0("Saving scaled env. suitability to: ", scaledSuitPath))

    return(globalMaxSuitability) # Return the global maximum suitability for scaling
}

#' Compute Environmental Suitability Using the Beta Function
#'
#' This function calculates environmental suitability scores based on the asymmetric beta function
#' described by Yan and Hunt (1999), which is a simplified version of the function originally proposed
#' by Yin et al. (1995). The suitability is 0 outside the defined environmental range
#' (between MinEnvVar and MaxEnvVar) and peaks at OptEnvVar.
#'
#' The formula used is:
#' \deqn{
#'   suitability = \left( \frac{V_{max} - V_{env}}{V_{max} - V_{opt}} \right)
#'                 \cdot \left( \frac{V_{env} - V_{min}}{V_{opt} - V_{min}} \right)^{\frac{V_{opt} - V_{min}}{V_{max} - V_{opt}}}
#' }
#'
#' where:
#' - \eqn{V_{env}} is the environmental value at a given time
#' - \eqn{V_{min}}, \eqn{V_{opt}}, and \eqn{V_{max}} are the minimum, optimum, and maximum values for suitability
#'
#' @param MinEnvVar Array minimum tolerated environmental values (no. species x no. env. variables).
#' @param MaxEnvVar Array maximum tolerated environmental values (no. species x no. env. variables).
#' @param OptEnvVar Array optimal environmental values (no. species x no. env. variables).
#' @param EnvVar A numeric array of actual environmental values (length x depth x height x env. variables).
#'
#' @return A numeric array of shape length x depth x height x no. species x no. env. variables,
#'         with suitability values in the range [0, 1].
#'
#' @references
#' Yan, Weikai, and L. A. Hunt (1999). An equation for modelling the temperature response of plants using only
#' the cardinal temperatures. *Annals of Botany*, 84(5), 607–614. \doi{10.1006/anbo.1999.0955}
#'
#' Yin, X., Kropff, M. J., McLaren, G., & Visperas, R. M. (1995). A nonlinear model for crop development
#' as a function of temperature. *Agricultural and Forest Meteorology*, 77(1-2), 1–16.
#'
#' @examples
#' # Simple example with arrays
#' MinEnvVar <- array(14, dim = c(100, 2)) # 100 species, 2 environmental variables
#' MaxEnvVar <- array(29, dim = c(100, 2))
#' OptEnvVar <- array(21, dim = c(100, 2))
#' EnvVar <- array(rnorm(50 * 50 * 60 * 2, mean = 21, sd = 12), dim = c(50, 50, 60, 2))
#' SuitabilityScore(MinEnvVar, MaxEnvVar, OptEnvVar, EnvVar)
SuitabilityScore <- function (MinEnvVar, MaxEnvVar, OptEnvVar, EnvVar) {
    # Dimensions
    spatial_dim <- dim(EnvVar)[1:3]  # [50, 50, 60]
    n_species <- dim(MinEnvVar)[1]  # 100
    n_vars <- dim(MinEnvVar)[2]     # 2

    # Expand EnvVar to [50, 50, 60, 100, 2]
    EnvVar_exp <- array(EnvVar, dim = c(dim(EnvVar), 1))
    EnvVar_exp <- array(EnvVar_exp, dim = c(dim(EnvVar), n_species))
    EnvVar_exp <- aperm(EnvVar_exp, c(1, 2, 3, 5, 4))

    # Expand Min/Opt/MaxEnvVar to [50, 50, 60, 100, 2]
    MinEnvVar_exp <- array(rep(MinEnvVar, each = prod(spatial_dim)),
                           dim = c(spatial_dim, n_species, n_vars))
    MaxEnvVar_exp <- array(rep(MaxEnvVar, each = prod(spatial_dim)),
                           dim = c(spatial_dim, n_species, n_vars))
    OptEnvVar_exp <- array(rep(OptEnvVar, each = prod(spatial_dim)),
                           dim = c(spatial_dim, n_species, n_vars))

    # Create zero array for output
    suitability <- array(0.0, dim = c(spatial_dim, n_species, n_vars))  # [50, 50, 60, 100, 2]

    # Valid mask: within bounds
    valid <- (EnvVar_exp > MinEnvVar_exp) & (EnvVar_exp < MaxEnvVar_exp)

    ValidMaxEnvVar_exp <- MaxEnvVar_exp[valid]
    ValidOptEnvVar_exp <- OptEnvVar_exp[valid]
    ValidMinEnvVar_exp <- MinEnvVar_exp[valid]
    ValidEnvVar <- EnvVar_exp[valid]

    # Pre-compute denominators
    MaxOptDiff <- ValidMaxEnvVar_exp - ValidOptEnvVar_exp
    OptMinDiff <- ValidOptEnvVar_exp - ValidMinEnvVar_exp

    # Compute suitability only for valid entries
    num <- (ValidMaxEnvVar_exp - ValidEnvVar) / MaxOptDiff
    denom <- (ValidEnvVar - ValidMinEnvVar_exp) / OptMinDiff
    expo  <- OptMinDiff / MaxOptDiff

    if (anyNA(valid)) {
      valid[is.na(valid)] <- TRUE # Make sure NA values are stored (i.e., they are no NAs in the mask)
    }
    suitability[valid] <- num * denom^expo

    return(suitability)  # shape: e.g. [50, 50, 60, 100, 2]
}


# Parabolic Optimum function
Parabol <- function(a, b, c, x) {

    if (length(a) > 1 & length(b) > 1 & length(c) > 1 & length(x) > 1) {
        n_species <- length(a)
        spatial_dim <- length(x)

        # Expand dimensions to match -> (length, depth, height, n_species)
        a_exp <- array(a, dim = c(1, 1, 1, n_species))
        a_exp <- array(a_exp, dim = c(spatial_dim, n_species))
        b_exp <- array(b, dim = c(1, 1, 1, n_species))
        b_exp <- array(b_exp, dim = c(spatial_dim, n_species))
        c_exp <- array(c, dim = c(1, 1, 1, n_species))
        c_exp <- array(c_exp, dim = c(spatial_dim, n_species))
        x_exp <- array(x, dim = c(spatial_dim, 1))
        x_exp <- array(x_exp, dim = c(spatial_dim, n_species))

        return((a_exp * x_exp^2) + (b_exp * x_exp) + c_exp)

    } else {
        return((a * x^2) + (b * x) + c)
    }
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
    MicrohabitatType <- config$MicrohabitatType  # Define which type of forest the microhabitat belongs to. 1: dynamic forest, 2: static forest, 3: uniform forest
    timeSteps <- config$timeSteps  # Model for timeSteps beginning at the time step given by the initial distribution
    InitialTimeStep <- config$InitialTimeStep
    # Choose species pools to use and number of replicates per species pool
    numSpeciesPools <- config$numSpeciesPools
    # Specifying light response function for growth
    LightResponseFct <- config$LightResponseFct
    numReplicate <- config$numReplicate
    Imax <- config$Imax  # Maximum light intensity

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
    if (LightResponseFct == "Yan and Hunt") {
        EnvScoreVars <- c("Hum", "Temp", "Wind", "Light")
    } else {
        EnvScoreVars <- c("Hum", "Temp", "Wind")
    }
    MhIdx <- Inds[paste0(EnvScoreVars, "NicheOpt")]

    # Create folder to save the model results
    dir.create(DirectoryOutput, recursive=TRUE)

    # Parallelize across species pools and species
    output <- foreach (numPool=seq(numSpeciesPools),
                       .export=c("ComputeSuitabilityUnscaled", "int_seq", "Parabol", "SuitabilityScore")) %dorng% {

        # First step: create probability matrices for each species
        # Load species pool
        SpeciesPoolFileName <- paste0("SpeciesPool", numPool, ".csv")
        SpeciesPool <- read.csv(file.path(DirectorySpeciesPools, SpeciesPoolFileName), sep=",", header=TRUE)

        # Create Save-Directory for each each replicate/initialDistribution
        DirectoryOutputSpeciesPool <- file.path(DirectoryOutput, paste0("ID_SpeciesP_", numPool, "_Rep_", numReplicate))
        dir.create(DirectoryOutputSpeciesPool, recursive=TRUE)

        globalMaxSuitability <- ComputeSuitabilityUnscaled(
          timeSteps,
          DirectoryMicrohabitat,
          SpeciesPool,
          InitialTimeStep,
          Imax,
          LightResponseFct,
          MhIdx,
          EnvScoreVars,
          DirectoryOutputSpeciesPool,
          dimPlot
        )
  }
}

main()