options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("utils.R")

library("foreach")
library("doParallel")
library("doRNG")

###############################################################################


fill_edge_na <- function(array3d, startX, startY, startZ, endX, endY, endZ) {
    dims <- dim(array3d)

    for (i in seq_len(dims[1])) {
        for (j in seq_len(dims[2])) {
            for (k in seq_len(dims[3])) {
                if (is.na(array3d[i, j, k])) {
                    # Clamp indices to valid range inside original wind matrix
                    ci <- min(max(i, startX), endX)
                    cj <- min(max(j, startY), endY)
                    ck <- min(max(k, startZ), endZ)

                    array3d[i, j, k] <- array3d[ci, cj, ck]
                }
            }
        }
    }
    return(array3d)
}


compute_prob_matrix_norm <- function(centralPoint, dimPlot, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool, WindSpeed) {
    # Erzeugen der Distanzmatrix mit allen Distanzen zum
    DistanceMatrix <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))

    for (i in seq_len(dimX)) {
        for (j in seq_len(dimY)) {
            for (k in seq_len(dimZ)) {
                x1 <- c(i, j, k)
                x2 <- c(centralPoint[1], centralPoint[2], centralPoint[3])
                DistanceMatrix[i, j, k] <- sqrt(sum((x1 - x2)^2))  # call to pdist() in the matlab script
            }
        }
    }

    # Modify the windspeed matrix to match the larger distance matrix
    WindSpeedFull <- array(NA, dim = c(dimX, dimY, dimZ))

    # Coordinates to insert original wind field in center
    startX <- floor(dimX / 2) - floor(dimPlot[1] / 2) + 1
    startY <- floor(dimY / 2) - floor(dimPlot[2] / 2) + 1
    startZ <- floor(dimZ / 2) - floor(dimPlot[3] / 2) + 1
    endX <- startX + dimPlot[1] - 1
    endY <- startY + dimPlot[2] - 1
    endZ <- startZ + dimPlot[3] - 1

    # Embed the original wind field into the center
    WindSpeedFull[startX:endX,
                  startY:endY,
                  startZ:endZ] <- WindSpeed

    # Erzeugen der Wahrschienlichkeitsmatrix anhand der Distanzmatrix und dem
    # artspezischien Wert bb aus der Epiphytenmatrix
    # Google translate: Generating the probability matrix based on the distance
    # matrix and the species value from the epiphyte matrix
    # negExp = @(distance,bb) exp(-distance.*bb); %Negative Exponential function
    ProbabilityMatrix <- array(rep(0, dimX * dimY * dimZ * NumberOfSpecies), dim=c(dimX, dimY, dimZ, NumberOfSpecies))
    ProbabilityMatrixNormalized <- array(rep(0, dimX * dimY * dimZ * NumberOfSpecies), dim=c(dimX, dimY, dimZ, NumberOfSpecies))

    for (i in seq_len(NumberOfSpecies)) {
        exponentE <- SpeciesPool$DispersalKernel[i]
        dispersalAsymmetry <- SpeciesPool$DispersalKernelAsymmetry[i]
        dispersalWindEffect <- SpeciesPool$DispersalKernelWindEffect[i]

        # Scale the dispersal kernel by wind speed (depending on species specific wind dispersal)
        WindExponentE <- exponentE / (1 + dispersalWindEffect * WindSpeedFull)
        # Replace NAs with former exponentE
        WindExponentE[is.na(WindExponentE)] <- exponentE

        ProbabilityMatrix[, , , i] <- exp(-DistanceMatrix * WindExponentE)  # call to negExp(DistanceMatrix(:,:,:),exponentE) in matlab

        print(paste0("exponentE for species ", i, ": ", mean(WindExponentE)))
        print(paste0("dispersalAsymmetry for species ", i, ": ", dispersalAsymmetry))
        print(paste0("Sum of ProbabilityMatrix for species ", i, ": ", sum(ProbabilityMatrix[, , , i])))

        # Apply dispersal asymmetry (probability to disperse downwards higher than upwards dispersal)
        # WARNING: The index i in the 4th dimension was missing from the Matlab script and so we were
        # getting "incorrect number of dimensions" in R. Added it here but need to ask if this is what
        # it was supposed to do.
        id3_1 <- int_seq(from=centralPoint[3], to=dimZ, by=1)
        id3_2 <- int_seq(from=1, to=centralPoint[3] - 1, by=1)
        ProbabilityMatrix[, , id3_1, i] <- ProbabilityMatrix[, , id3_1, i] * ((1 - dispersalAsymmetry) / 0.5)
        ProbabilityMatrix[, , id3_2, i] <- ProbabilityMatrix[, , id3_2, i] * (dispersalAsymmetry / 0.5)

        ProbabilityMatrixNormalized[, , , i] <- ProbabilityMatrix[, , , i] / sum(ProbabilityMatrix[, , , i])  # sum(sum(sum(ProbabilityMatrix(:,:,:,i)))) in matlab
    }

    return(ProbabilityMatrixNormalized)
}


ComputeSuitabilityUnscaled <- function (timeSteps,
                                   DirectoryMicrohabitat,
                                   SpeciesPool,
                                   InitialTimeStep,
                                   Imax,
                                   LightResponseFct,
                                   Inds,
                                   DirectoryOutput,
                                   Overwrite = FALSE) {
    NSpecies <- nrow(SpeciesPool)
    globalMaxSuitability <- rep(-Inf, NSpecies)

    for (t in seq_len(timeSteps)) {

        envSuitPath <- file.path(DirectoryOutput, paste0("EnvSuitability_t", InitialTimeStep + t - 1, ".rds"))

        if (!Overwrite & file.exists(envSuitPath)) { # Skip iteration if EnvSUitability already exists
            next
        }

        MHPath <- file.path(DirectoryMicrohabitat,
                            paste0("MicrohabitatMatrix", InitialTimeStep + t - 1, ".rds"))
        Microhabitat <- readRDS(MHPath)
        LightIdx <- Inds["LightNicheOpt"]
        Microhabitat[, , , LightIdx] <- Microhabitat[, , , LightIdx] * Imax

        # Check which light response function to use and set variables for suitability calculation
        if (LightResponseFct == "Yan and Hunt") {
            EnvScoreVars <- c("Hum", "Temp", "Wind", "Light")
        } else {
            EnvScoreVars <- c("Hum", "Temp", "Wind")
        }

        MhIdx <- Inds[paste0(EnvScoreVars, "NicheOpt")]

        # Extract environmental variables from the species pool
        MinEnvVar <- as.matrix(SpeciesPool[paste0("Min", EnvScoreVars)])
        MaxEnvVar <- as.matrix(SpeciesPool[paste0("Max", EnvScoreVars)])
        OptEnvVar <- as.matrix(SpeciesPool[paste0("Optimum", EnvScoreVars)])

        # -- Compute the score for each environmental variable
        spatialDim <- dim(Microhabitat)
        # Initialize the array for environmental suitability scores with zeros
        EnvSuitabilityVars <- array(
          rep(0, spatialDim[1] * spatialDim[2] * spatialDim[3] * NSpecies * 4),
          dim=c(spatialDim[1], spatialDim[2], spatialDim[3], NSpecies, 4))
        # Compute score for either Hum, Temp, Wind, and Light or only Hum, Temp, and Wind
        EnvSuitabilityVars[, , , , 1:length(EnvScoreVars)] <- SuitabilityScore(
          MinEnvVar, MaxEnvVar, OptEnvVar, Microhabitat[ , , , MhIdx])

        if (LightResponseFct == "Parabolic") {
            EnvSuitabilityVars[, , , , 4] <- Parabol(
              SpeciesPool$LightResponseA, SpeciesPool$LightResponseB,
              SpeciesPool$LightResponseC, Microhabitat[ , , , 3])
        }

        # Combine the suitability probabilities for all environmental variables
        EnvSuitability <- apply(EnvSuitabilityVars, c(1, 2, 3, 4), prod)

        # Get the maximum suitability for this time step for later scaling
        maxThisStep <- apply(EnvSuitability, c(4), max, na.rm = TRUE)
        isNewMax <- maxThisStep > globalMaxSuitability
        globalMaxSuitability[isNewMax] <- maxThisStep[isNewMax]

        # Save to disk
        saveRDS(EnvSuitability, envSuitPath)
    }

    return(globalMaxSuitability) # Return the global maximum suitability for scaling
}


# Functions used in the model

# Bertalanffy Growth
GrowthRate <- function(MaxMass, Mass, K) {
    return(K * (MaxMass - Mass))
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

    # Compute suitability only for valid entries
    num   <- (MaxEnvVar_exp[valid] - EnvVar_exp[valid]) / (MaxEnvVar_exp[valid] - OptEnvVar_exp[valid])
    denom <- (EnvVar_exp[valid] - MinEnvVar_exp[valid]) / (OptEnvVar_exp[valid] - MinEnvVar_exp[valid])
    expo  <- (OptEnvVar_exp[valid] - MinEnvVar_exp[valid]) / (MaxEnvVar_exp[valid] - OptEnvVar_exp[valid])

    valid[is.na(valid)] <- TRUE # Make sure NA values are stored (i.e., they are no NAs in the mask)
    suitability[valid] <- num * denom^expo

    return(suitability)  # shape: [50, 50, 60, 100, 2]
}


dispersal <- function(NumberOfSpecies,
                      E,
                      Microhabitat,
                      SurfaceBiomassScaling,
                      dimPlot,
                      centralPoint,
                      InterceptRecruitment,
                      SlopeRecruitment,
                      ProbabilityMatrixNormalized,
                      SpeciesPool,
                      MaxIndividualID,
                      Inds) {
    # Store number of individuals at beginning of time step
    IntialNumberIndividuals <- array(rep(0, NumberOfSpecies))
    for (g in seq_len(NumberOfSpecies)) {
        # Count indices where SpeciesID is g and Status is 1
        IntialNumberIndividuals[g] <- length(which(E$SpeciesID == g & E$Status == 1))
    }
    IntialNumberIndividualsTotal <- length(which(E$Status == 1))
    InitialNumberSpecies <- length(unique(E$SpeciesID[E$Status == 1]))
    NumberRecruitsPerSpecies <- array(rep(0, NumberOfSpecies))

    # Calculate free surface area per voxel
    AvailableSurfaceArea <- Microhabitat[, , , Inds["TotalSurfaceAreaOpt"]]
    for (i in seq_len(nrow(E))) {
        SurfaceAreaNeededInVoxel <- E$Mass[i]^(2/3) / SurfaceBiomassScaling
        AvailableSurfaceArea[E$X[i], E$Y[i], E$Z[i]] <- max(0, AvailableSurfaceArea[E$X[i], E$Y[i], E$Z[i]] - SurfaceAreaNeededInVoxel)
    }

    # Check if there are species left (~isempty(E) in matlab)
    if (nrow(E) > 0) {
        unique_species <- unique(E$SpeciesID)  # list with species IDs of all present species

        # Initialize potential recruitment dataframe
        PotentialRecruitment <- data.frame(matrix(0, nrow=length(unique_species), ncol=2))
        colnames(PotentialRecruitment) <- c("index", "potential_recruit")

        # loop over all species
        for (i in seq_len(length(unique_species))) {
            # Generate initially empty matrix to store the probabilities for recruitment
            ProbabilityMatrixPerSpecies <- array(rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]), dim=c(dimPlot[1], dimPlot[2], dimPlot[3]))

            # Matrix containing all mature individuals of one species
            MatureIndividulsPerSpecies <- E[E$SpeciesID == unique_species[i] & E$Mass >= E$MassAtMaturity, ]

            # ~isempty(MatureIndividulsPerSpecies) in matlab
            if (nrow(MatureIndividulsPerSpecies) > 0) {

                # Probability matrix for each species: Depending on the position of each mature individual,
                # the total probability for the species is calculated.
                # The second part of the equation accounts for the actual size of the individual
                # in relation to the maximum size for which the recruitment per individual is defined
                for (j in seq_len(nrow(MatureIndividulsPerSpecies))) {
                    idx1 <- seq(from=centralPoint[1] - MatureIndividulsPerSpecies$X[j] + 1, to=centralPoint[1] - MatureIndividulsPerSpecies$X[j] + dimPlot[1], by=1)
                    idx2 <- seq(from=centralPoint[2] - MatureIndividulsPerSpecies$Y[j] + 1, to=centralPoint[2] - MatureIndividulsPerSpecies$Y[j] + dimPlot[2], by=1)
                    idx3 <- seq(from=centralPoint[3] - MatureIndividulsPerSpecies$Z[j] + 1, to=centralPoint[3] - MatureIndividulsPerSpecies$Z[j] + dimPlot[3], by=1)
                    idx4 <- MatureIndividulsPerSpecies$SpeciesID[j]

                    factor1 <- (InterceptRecruitment + (SlopeRecruitment * MatureIndividulsPerSpecies$Mass[j])) * MatureIndividulsPerSpecies$RecruitmentInvestmentRel[j]
                    factor2 <- (MatureIndividulsPerSpecies$Mass[j] - MatureIndividulsPerSpecies$MassAtMaturity[j]) / (MatureIndividulsPerSpecies$MaximumMass[j] - MatureIndividulsPerSpecies$MassAtMaturity[j])
                    factor3 <- 1 + (MatureIndividulsPerSpecies$RecruitmentInc[j] * factor2)

                    ProbabilityMatrixPerSpecies <- ProbabilityMatrixPerSpecies + ProbabilityMatrixNormalized[idx1, idx2, idx3, idx4] * factor1 * factor3
                }

                # Store potential normalized number of recruits. We will use this to populate SummaryMatrixSpecies later
                PotentialRecruitment$index[i] <- i
                PotentialRecruitment$potential_recruit[i] <- sum(ProbabilityMatrixPerSpecies)  # potential recruitment / sum(sum(sum(ProbabilityMatrixPerSpecies))) in matlab

                # Matix containing all voxel for which the light requirements are fulfilled
                # We use the first row from MatureIndividulsPerSpecies. Since its elements have the same SpeciesID
                # then the MinLight and MaxLight is the same for all rows.
                IdxLight <- Inds["LightNicheOpt"]
                IdxHum <- Inds["HumNicheOpt"]
                IdxTemp <- Inds["TempNicheOpt"]
                IdxWind <- Inds["WindNicheOpt"]

                LightSuitable <- ((Microhabitat[, , , IdxLight] >= MatureIndividulsPerSpecies$MinLight[1]) &
                                  (Microhabitat[, , , IdxLight] <= MatureIndividulsPerSpecies$MaxLight[1]))
                HumSuitable <- ((Microhabitat[, , , IdxHum] >= MatureIndividulsPerSpecies$MinHum[1]) &
                                (Microhabitat[, , , IdxHum] <= MatureIndividulsPerSpecies$MaxHum[1]))
                TempSuitable <- ((Microhabitat[, , , IdxTemp] >= MatureIndividulsPerSpecies$MinTemp[1]) &
                                 (Microhabitat[, , , IdxTemp] <= MatureIndividulsPerSpecies$MaxTemp[1]))
                WindSuitable <- ((Microhabitat[, , , IdxWind] >= MatureIndividulsPerSpecies$MinWind[1]) &
                                 (Microhabitat[, , , IdxWind] <= MatureIndividulsPerSpecies$MaxWind[1]))
                pot_habitat <- LightSuitable & HumSuitable #& TempSuitable & WindSuitable

                print(paste("Potential habitat no. voxels:", sum(pot_habitat, na.rm=TRUE)))

                # Final probabiliy matrix for new recruits
                probability_recruits <- ProbabilityMatrixPerSpecies * AvailableSurfaceArea * pot_habitat

                # Calculate number of recuits based on final probability matrix
                Recruits <- array(rpois(length(probability_recruits), probability_recruits), dim=dim(probability_recruits))  # poissrnd(probability_recruits) in matlab

                # Add new recruits to epiphyte matrix
                num_recruits <- sum(Recruits, na.rm=TRUE)  # sum(sum(sum(Recruits))) in matlab
                NumberRecruitsPerSpecies[unique_species[i]] <- num_recruits

                if (num_recruits > 0) {
                    ids <- arrayInd(which(Recruits > 0), dim(Recruits))
                    xInd <- ids[, 1]
                    yInd <- ids[, 2]
                    zInd <- ids[, 3]

                    while (num_recruits > length(xInd)) {
                        tmp_ids <- arrayInd(which(Recruits > 0), dim(Recruits))
                        Recruits[tmp_ids] = Recruits[tmp_ids] - 1

                        tmp_ids <- arrayInd(which(Recruits > 0), dim(Recruits))
                        xInd <- append(xInd, tmp_ids[, 1])
                        yInd <- append(yInd, tmp_ids[, 2])
                        zInd <- append(zInd, tmp_ids[, 3])
                    }
                    vec_recruits <- seq(from=nrow(E) + 1, to=nrow(E) + length(xInd), by=1)

                    # Copy species information to Epiphyte matrix
                    E[vec_recruits, names(SpeciesPool)] <- SpeciesPool[unique_species[i], ]
                    E$X[vec_recruits] <- xInd
                    E$Y[vec_recruits] <- yInd
                    E$Z[vec_recruits] <- zInd
                    E$Mass[vec_recruits] <- 0  # Initial size
                    E$Status[vec_recruits] <- 1  # status 1:alive
                    E$IndividualID[vec_recruits] <- seq(from=MaxIndividualID + 1, to=MaxIndividualID + length(xInd), by=1)  # individual ID

                    E[is.na(E)] <- 0  # convert all NA to 0 so that the R script matches the Matlab

                    MaxIndividualID <- MaxIndividualID + length(xInd)
                }
            }
        }
    } else {
        # Initialize empty potential recruitment dataframe
        PotentialRecruitment <- data.frame(matrix(0, nrow=0, ncol=2))
        colnames(PotentialRecruitment) <- c("index", "potential_recruit")
    }

    disp_items <- list("IntialNumberIndividuals"=IntialNumberIndividuals,
                       "NumberRecruitsPerSpecies"=NumberRecruitsPerSpecies,
                       "InitialNumberSpecies"=InitialNumberSpecies,
                       "IntialNumberIndividualsTotal"=IntialNumberIndividualsTotal,
                       "E"=E,
                       "PotentialRecruitment"=PotentialRecruitment,
                       "MaxIndividualID"=MaxIndividualID)
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
    if (nTasks != "") {
        numCores <- strtoi(nTasks)
    }
    else {
        numCores <- detectCores() - 1
    }

    registerDoParallel(numCores)

    # Parse input configuration file
    config <- parse_config()

    ###############################################################################
    # Parameters that need to be specified/checked before running this script

    # Input directories
    DirectoryMicrohabitat <- config$DirectoryMicrohabitat
    DirectorySpeciesPools <- config$DirectorySpeciesPools
    DirectoryModelMain <- config$DirectoryModelMain

    # Output directory
    DirectoryModelResults <- config$DirectoryModelResults

    MicrohabitatType <- config$MicrohabitatType  # Define which type of forest the microhabitat belongs to. 1: dynamic forest, 2: static forest, 3: uniform forest

    # Model parameters
    timeSteps <- config$timeSteps  # Model for timeSteps beginning at the time step given by the initial distribution

    # Density of individuals per ha at which to stop the simulationof the community and
    # move to the next replicate (to prevent exploding communities)
    StopCriterionHa <- config$StopCriterionHa  # Individuals per ha

    # Choose species pools to use and number of replicates per species pool
    numSpeciesPools <- config$numSpeciesPools  # Start and end number of  species pools (if the species pools do not exist, they are automatically skipped)
    replicatePerSpeciesPool <- config$replicatePerSpeciesPool  # Number of replicates per species pool  (if the replicates do not exist, they are automatically skipped)

    SurfaceBiomassScaling <- config$SurfaceBiomassScaling  # cm^2 per m^2
    Imax <- config$Imax  # maximum light above canopy

    # Competition Methods; defines which individuals are removed in voxels which
    # are entirely filled. 1:size (small individuals are outcompetet by larger ones); 2:random competition
    CompetitionMethod <- config$CompetitionMethod

    # Specifying light response function for growth
    LightResponseFct <- config$LightResponseFct

    # Mortality method (complete random or scaling with mass according to metabolic theory);
    MortalityMethod <- config$MortalityMethod  # 0: random mortality; 1: scaling with mass to the exponent -1/4
    MortRateRandom <- config$MortRateRandom
    MortRateMass <- config$MortRateMass
    MortRateMassScaling <- config$MortRateMassScaling  # widely used scaling fator

    InitialTimeStep <- config$InitialTimeStep  # Time step for which the Initial distribution is generated in A3

    # RNG seed
    random_state_file <- config$RandomState
    if (!is.null(random_state_file) && file.exists(random_state_file)) {
        writeLines("Loading previous random number generator state...")
        restore_rng(random_state_file)
    } else {
        # If we provide an integer, use it to set the seed
        # otherwise it will be NULL and therefore a random
        # seed will be created.
        seed <- config$seed
        set.seed(seed, kind="Mersenne-Twister")
    }
    
    # Get microhabitat matrix variables - dynamic handling of selected variables
    MicrohabitatVariableFlags <- config$MicrohabitatVariableFlags
    MhVarNames <- c("TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt", "AverageWeightedAngles",
                    "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
    # Only keep active options
    ActiveOpts <- MhVarNames[as.logical(MicrohabitatVariableFlags)]
    # Assign indices
    Inds <- setNames(seq_along(ActiveOpts), ActiveOpts)
    
    ###############################################################################

    # Create folder to save the model results
    dir.create(DirectoryModelResults, recursive=TRUE)

    # Save current random state.
    # We need to store this after setting the seed but before calling any
    # functions that generate random numbers.
    writeLines("Storing the random number generator state...")
    save_rng(file.path(DirectoryModelResults, "random_state_seed.RData"))

    # Load plot dimensions
    dimPlot <- readRDS(file.path(DirectoryMicrohabitat, "dimPlot.rds"))

    # Set StopCriterion for this simulation
    StopCriterion <- 0.0001 * StopCriterionHa * dimPlot[1] * dimPlot[2]

    # Load TraitRanges (ranges used to create the species pool)
    FileTraitRanges <- file.path(DirectorySpeciesPools, "TraitRanges.csv")
    TraitRanges <- read.table(FileTraitRanges, sep=",", header=FALSE)

    SlopeRecruitment <- TraitRanges[1, 1]
    InterceptRecruitment <- TraitRanges[2, 1]


    # Following the columns in the species matrix refering to a trait or
    # variable. This is handy if the epiphyte matrix changes
    ColSSpeciesID <- 1
    ColSNumberIndividualsBeginning <- 2
    ColSNumberIndividualsEnd <- 3
    ColSNumberMatureIndividuals <- 4
    ColSNumberRecruits <- 5
    ColSNumberRecruitsPotential <- 6
    ColSNumberMortalityBranchFall <- 7
    ColSNumberMortalityLight <- 8
    ColSNumberMortalityCompetition <- 9
    ColSNumberMortalityNatural <- 10
    ColSNumberPopulationGrowthRate <- 11
    ColSNumberPopulationGrowthRateLog <- 12
    ColSNumberBirthRate <- 13
    ColSNumberDeathRate <- 14
    ColSAverageSize <- 15
    ColSAverageAge <- 16
    ColSMinLight <- 17
    ColSMaxLight <- 18
    ColSMeanLight <- 19
    # ColSMinHeight <- 20
    # ColSMaxHeight <- 21
    # ColSMeanHeight <- 22
    # Climate mortality columns
    ColSNumberMortalityHum <- 20
    ColSNumberMortalityTemp <- 21
    ColSNumberMortalityWind <- 22
    ColSMinHum <- 23
    ColSMaxHum <- 24
    ColSMeanHum <- 25
    ColSMinTemp <- 26
    ColSMaxTemp <- 27
    ColSMeanTemp <- 28
    ColSMinWind <- 29
    ColSMaxWind <- 30
    ColSMeanWind <- 31

    TotalColsSpeciesMatrix <- ColSMeanWind

    # Headers of matrix
    SummaryMatrixSpeciesHeaders <- c("TimeStep", "SpeciesID", "NumberIndividualsBeginning", "NumberIndividualsEnd",
        "NumberMatureIndividuals", "NumberRecruits", "NumberRecruitsPotential", "NumberMortalityBranchFall",
        "NumberMortalityLight", "NumberMortalityCompetition", "NumberMortalityNatural", "PopulationGrowthRate",
        "PopulationGrowthRateLog", "BirthRate", "DeathRate", "AverageMass", "AverageAge", "MinLight", "MaxLight",
        "MeanLight", #"MinHeight", "MaxHeight", "MeanHeight"
        "NumberMortalityHum", "NumberMortalityTemp", "NumberMortalityWinds", "MinHum", "MaxHum", "MeanHum",
        "MinTemp", "MaxTemp", "MeanTemp", "MinWind", "MaxWind", "MeanWind")

    if (length(SummaryMatrixSpeciesHeaders) != (TotalColsSpeciesMatrix + 1)) {
        stop("Headers of species matrix do not match with number of columns")
    }

    # Headers of matrix
    SummaryMatrixCommunityHeaders <- c("timeStep", "NumberSpeciesBeginning", "NumberSpeciesEnd",
        "NumberIndividualsBeginning", "NumberIndividualsEnd", "Recruits", "MortalityBranchFall",
        "MortalityLight", "MortalityCompetition", "MortalityNatural", "MortalityHum", "MortalityTemp",
        "MortalityWind", "BranchSurfaceIndex", "EpiphyteFilling")

    ###############################################################################
    # Main loop for the community model for each species pool and for each replicate
    pairs <- create_pairs(numSpeciesPools, replicatePerSpeciesPool)

    # Each loop employs a different random number generator (RNG) stream, resulting in distinct,
    # statistically independent random sequences. These sequences are reproducible across multiple
    # runs, provided the master seed and other input parameters remain unchanged. Importantly, the
    # number of cores does not influence the sequences, only the order in which the loops are
    # processed. Consequently, parallel and serial execution will yield identical results.
    # Internally, the foreach package employs the L'Ecuyer-CMRG RNG algorithm for reliable random
    # number generation, ensuring reproducible results even in parallel computing environments.
    output <- foreach (pair_idx=seq_len(nrow(pairs)),
                       .export=c("ComputeSuitabilityUnscaled", "compute_prob_matrix_norm", "int_seq", "dispersal", "GrowthRate", "Parabol")) %dorng% {
        numPool <- pairs$numPool[pair_idx]
        r <- pairs$r[pair_idx]

        # Check if a initial distribution for the species pool exists. If not, move on to the next species pool
        FileNameInitalDistribution <- file.path(DirectoryModelMain, paste0("ID_SpeciesP_", numPool, "_Rep_", r, ".csv"))
        if (!file.exists(FileNameInitalDistribution)) {
            print(paste0("Initial distribution file ", FileNameInitalDistribution,
                         " does not exist. Skipping species pool ", numPool, ", replicate ",
                         r, "."))
            return(NULL)
        }

        # First step: create probability matrices for each species
        # Load species pool
        SpeciesPoolFileName <- paste0("SpeciesPool", numPool, ".csv")
        SpeciesPool <- read.csv(file.path(DirectorySpeciesPools, SpeciesPoolFileName), sep=",", header=TRUE)
        NumberOfSpecies <- nrow(SpeciesPool)  # number of species per 25X25m plot

        ###########################################################################
        # Erzeugen der Distanzmatrix und der Wahrscheinlichkeitsmatrix für jede Art
        # Dimensionen der Dispersal matrix
        dimX <- dimPlot[1] * 2 + 1
        dimY <- dimPlot[2] * 2 + 1
        dimZ <- dimPlot[3] * 2 + 1

        centralPoint <- c(floor(dimX/2) + 1, floor(dimY/2) + 1, floor(dimZ/2) + 1)

        # Create Save-Directory for each each replicate/initialDistribution
        DirectoryModelResultsRun <- file.path(DirectoryModelResults, paste0("ID_SpeciesP_", numPool, "_Rep_", r))
        dir.create(DirectoryModelResultsRun, recursive=TRUE)

        if (MicrohabitatType == 1) {  # Dynamic forest
            # Filter species pool to only include initialized species

            # Precompute unscaled env. suitabilities and get the global maximum suitability for scaling
            SpeciesMaxSuitability <- ComputeSuitabilityUnscaled(
              timeSteps, DirectoryMicrohabitat, SpeciesPool, InitialTimeStep, Imax,
              LightResponseFct, Inds, DirectoryModelResultsRun, Overwrite = FALSE)
        }

        # Load initial epiphyte distribution
        E <- read.csv(FileNameInitalDistribution, sep=",", header=TRUE)  # E for epiphytes

        # Add column to E for additional information
        E[, c("TotalSurfaceInVoxel", "LightInVoxel", "HumInVoxel", "TempInVoxel", "WindInVoxel", "SurfaceLossInVoxel")] <- 0

        MaxIndividualID <- nrow(E)  # to trace individual IDs

        # Initialize Matrix where community parameters are save
        SummaryMatrixCommunity <- data.frame(matrix(0.0, nrow=timeSteps, ncol=length(SummaryMatrixCommunityHeaders)))
        colnames(SummaryMatrixCommunity) <- SummaryMatrixCommunityHeaders

        # Load microhabitat matrix if a uniform or static forest is simulated (only needs to be loaded once an not envery timestep)
        if (MicrohabitatType == 2 || MicrohabitatType == 3) {
            Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, "MicrohabitatMatrix1.rds"))
            Microhabitat[, , , Inds["LightNicheOpt"]] <- Microhabitat[, , , Inds["LightNicheOpt"]] * Imax  # In the microhabitat matrix, the realtive light extinction is stored: convert to light values in ?mol*m-2*s-1

            d1 <- dim(Microhabitat)[1]
            d2 <- dim(Microhabitat)[2]
            d3 <- dim(Microhabitat)[3]
            pot_habitat <- array(rep(0, d1 * d2 * d3), dim=c(d1, d2, d3))
        }

        # Initialize matrices where the aggregated information on species level are saved
        SummaryMatrixSpeciesSave <- array(rep(0, (timeSteps*NumberOfSpecies) * (TotalColsSpeciesMatrix + 1)), dim=c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix + 1))

        # Initialize Matrix where speceies parameters are save
        SummaryMatrixSpecies <- array(rep(0, (timeSteps*NumberOfSpecies) * TotalColsSpeciesMatrix), dim=c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix))

        for (t in seq_len(timeSteps)) {

            # Check if the stop criterion is met
            if (length(which(E$Status == 1)) > StopCriterion) {
                break
            }

            # Load microhabitat matrix for specific timeStep if dynamic forest is simulated
            if (MicrohabitatType == 1) {
                Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, paste0("MicrohabitatMatrix", InitialTimeStep + t - 1, ".rds")))
                Microhabitat[, , , Inds["LightNicheOpt"]] <- Microhabitat[, , , Inds["LightNicheOpt"]] * Imax  # In the microhabitat matrix, the realtive light extinction is stored: convert to light values in ?mol*m-2*s-1

                # Compute scaled suitability scores for the microhabitat
                file <- file.path(DirectoryModelResultsRun, paste0("EnvSuitability_t", t, ".rds"))
                if (file.exists(file)) {
                    EnvSuitability <- readRDS(file)
                    # Scale the suitability scores by the global specoes maximum suitability
                    SpeciesMaxSuitability_ext <- array(rep(SpeciesMaxSuitability, each = prod(dimPlot)),
                                                      dim = c(dimPlot, NumberOfSpecies))
                    ScaledEnvSuitability <- EnvSuitability / SpeciesMaxSuitability_ext
                } else {
                    stop("EnvSuitability file does not exist for time step ", t)
                }
            }

            ###############################################################################
            # 1. Dispersal

            # Create probability matrix for each species
            ProbabilityMatrixNormalized <- compute_prob_matrix_norm(
              centralPoint, dimPlot, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool,
              Microhabitat[,,,Inds["WindNicheOpt"]]
            )

            # Generate dispersal/recruitment matrix
            disp_items <- dispersal(
                NumberOfSpecies,
                E,
                Microhabitat,
                SurfaceBiomassScaling,
                dimPlot,
                centralPoint,
                InterceptRecruitment,
                SlopeRecruitment,
                ProbabilityMatrixNormalized,
                SpeciesPool,
                MaxIndividualID,
                Inds
            )

            # Out only.
            # Created in dispersal() and used later in the script.
            IntialNumberIndividuals <- disp_items$IntialNumberIndividuals
            NumberRecruitsPerSpecies <- disp_items$NumberRecruitsPerSpecies
            InitialNumberSpecies <- disp_items$InitialNumberSpecies
            IntialNumberIndividualsTotal <- disp_items$IntialNumberIndividualsTotal
            PotentialRecruitment <- disp_items$PotentialRecruitment

            # Inout.
            # Created outside dispersal(), modified in dispersal and used later too.
            E <- disp_items$E
            MaxIndividualID <- disp_items$MaxIndividualID  # This is only modified in dispersal()

            # Store potential normalized number of recruits in SummaryMatrixSpecies
            for (ii in seq_len(nrow(PotentialRecruitment))) {
                kk <- PotentialRecruitment$index[ii]
                if (kk != 0) {
                    SummaryMatrixSpecies[((kk-1) * timeSteps) + t, ColSNumberRecruitsPotential] <- PotentialRecruitment$potential_recruit[ii]
                }
            }

            NumberRecruits <- length(which(E$Status == 1)) - IntialNumberIndividualsTotal

            # Unclear what this line in the Matlab script is supposed to do.
            # From what I understand, the first column in E ("SpeciesID") takes non-zero values
            # only, so I think that E(:,1)==0 will always be empty.
            # E(E(:,1)==0,:)=[]; %in rare case, some individuals with only zeros are creates, which is wrong. This is to prevent the script to stop.

            ###############################################################################
            # Growth
            for (i in seq_len(nrow(E))) {
                # maybe it is faster if I do not use the if statement => speed testing
                if (E$Status[i] == 1) {
                    GrowthR <- GrowthRate(E$MaximumMass[i], E$Mass[i], E$GrowthRate[i])
                    SuitabilityIndividual <- ScaledEnvSuitability[E$X[i], E$Y[i], E$Z[i], E$SpeciesID[i]]

                    # TODO temprary remove later!
                    if (is.na(SuitabilityIndividual)) {
                        warning("SuitabilityIndividual is NA for individual ",
                                E$IndividualID[i],
                                ". Using species mean suitability for given height instead.")
                        SuitabilityIndividual <- mean(ScaledEnvSuitability[,,E$Z[i], E$SpeciesID[i]], na.rm = TRUE)
                        if (is.na(SuitabilityIndividual)) {
                            SuitabilityIndividual <- 0  # If still NA, set to 0
                        }
                    }

                    E$Mass[i] <- E$Mass[i] + max(0, GrowthR * SuitabilityIndividual)
                }

                # Add information about the voxel to the epiphyte matrix
                E$SurfaceAreaOccupied[i] <- (E$Mass[i]^(2/3)) / SurfaceBiomassScaling
                E$TotalSurfaceInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["TotalSurfaceAreaOpt"]]  # Total surface in voxel
                E$SurfaceLossInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["SurfaceAreaLossOpt"]]  # Percentage surface loss in this year
                E$LightInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["LightNicheOpt"]]  # Light conditions in voxel
                E$HumInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["HumNicheOpt"]]  # Humidity in voxel
                E$TempInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["TempNicheOpt"]]  # Temperature in voxel
                E$WindInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["WindNicheOpt"]]  # Wind in voxel
            }
            ###############################################################################

            ###############################################################################
            # Mortality
            for (i in seq_len(nrow(E))) {
                if (E$Status[i] == 1) {

                    # The following comparison would fail without the is.nan check,
                    # because Microhabitat contains NaNs in some entries and
                    # in R a comparison with a NaN returns NA, not a boolean.
                    # Note: We call runif repeatedly intentionally. See Issue #16 on Github
                    if (!is.nan(Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["SurfaceAreaLossOpt"]]) &&
                      runif(1, min=0, max=1) < Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["SurfaceAreaLossOpt"]]) {  # Mortality due to branch fall
                        E$Status[i] <- 3
                    } else if (E$LightInVoxel[i] < E$MinLight[i] | E$LightInVoxel[i] > E$MaxLight[i]) {  # Mortality due to changing light conditions
                        E$Status[i] <- 4
                    } else if (MortalityMethod == 0 && runif(1, min=0, max=1) < MortRateRandom) {  # Natural mortality rate
                        E$Status[i] <- 5
                    } else if (MortalityMethod == 1 && runif(1, min=0, max=1) < (MortRateMass * (E$Mass[i]^MortRateMassScaling))) {
                        E$Status[i] <- 5
                    } else if (!is.na(E$HumInVoxel[i]) && (E$HumInVoxel[i] < E$MinHum[i] | E$HumInVoxel[i] > E$MaxHum[i])) {
                        E$Status[i] <- 6
                    } else if (!is.na(E$TempInVoxel[i]) && (E$TempInVoxel[i] < E$MinTemp[i] | E$TempInVoxel[i] > E$MaxTemp[i])) {
                        E$Status[i] <- 7
                    } else if (!is.na(E$WindInVoxel[i]) && E$WindInVoxel[i] > E$MaxWind[i]) {
                        E$Status[i] <- 8
                    }
                }
            }

            # Mortality due to competition for space

            # Calculate total surface area occupied by epiphytes per voxel
            TotalSurfaceArePerVoxelOccupied <- array(rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]), dim=c(dimPlot[1], dimPlot[2], dimPlot[3]))
            for (w in seq_len(nrow(E))) {
                if (E$Status[w] == 1) {
                    TotalSurfaceArePerVoxelOccupied[E$X[w], E$Y[w], E$Z[w]] <- TotalSurfaceArePerVoxelOccupied[E$X[w], E$Y[w], E$Z[w]] + E$SurfaceAreaOccupied[w]
                }
            }

            # Indices of voxel where total area of epiphytes exeeds the available surface area
            ind_tmp <- arrayInd(which(TotalSurfaceArePerVoxelOccupied >
                                        Microhabitat[, , , Inds["TotalSurfaceAreaOpt"]]),
                                dim(TotalSurfaceArePerVoxelOccupied))
            IndX <- ind_tmp[, 1]
            IndY <- ind_tmp[, 2]
            IndZ <- ind_tmp[, 3]

            for (i in seq_len(length(IndX))) {
                # Get all epis in voxel
                EpisInVoxel <- E[E$X == IndX[i] & E$Y == IndY[i] & E$Z == IndZ[i] & E$Status == 1, ]

                # Sort them by size (CompetitionMethod=1) or randomly (CompetitionMethod=2)
                if (CompetitionMethod == 1) {
                    EpisInVoxel <- EpisInVoxel[order(EpisInVoxel$SurfaceAreaOccupied, decreasing=TRUE), ]
                } else if (CompetitionMethod == 2) {
                    EpisInVoxel <- EpisInVoxel[sample(seq_len(nrow(EpisInVoxel))), ]
                }

                CumulativeSumOfSurface <- cumsum(EpisInVoxel$SurfaceAreaOccupied)
                CumulativeSumOfSurfaceSum <- length(which(CumulativeSumOfSurface <= Microhabitat[IndX[i], IndY[i], IndZ[i], Inds["TotalSurfaceAreaOpt"]]))

                if (CumulativeSumOfSurfaceSum < nrow(EpisInVoxel)) {
                    E[is.element(E$IndividualID, EpisInVoxel[int_seq(CumulativeSumOfSurfaceSum + 1, nrow(EpisInVoxel)), "IndividualID"]), "Status"] <- 2
                }
            }
            ###############################################################################

            # Increase age
            E$Age <- E$Age + 1

            # Save number of mortality event
            MortalityCompetition <- length(which(E$Status == 2))
            MortalityBranchFall <- length(which(E$Status == 3))
            MortalityLight <- length(which(E$Status == 4))
            MortalityNatural <- length(which(E$Status == 5))
            MortalityHum <- length(which(E$Status == 6))
            MortalityTemp <- length(which(E$Status == 7))
            MortalityWind <- length(which(E$Status == 8))

            ###############################################################################

            # Store information in SummaryMatrixSpecies (summary over time for each species
            for (numSpecies in seq_len(NumberOfSpecies)) {

                # Define the row index once
                rowIndex <- ((numSpecies - 1) * timeSteps) + t

                SummaryMatrixSpecies[rowIndex, ColSSpeciesID] <- numSpecies
                SummaryMatrixSpecies[rowIndex, ColSNumberIndividualsBeginning] <- IntialNumberIndividuals[numSpecies]
                SummaryMatrixSpecies[rowIndex, ColSNumberIndividualsEnd] <- sum(E$Status == 1 & E$SpeciesID == numSpecies, na.rm=TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMatureIndividuals] <- sum(E$Status == 1 & E$SpeciesID == numSpecies & E$Mass >= E$MassAtMaturity, na.rm=TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberRecruits] <- NumberRecruitsPerSpecies[numSpecies]
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityBranchFall] <- sum(E$Status == 3 & E$SpeciesID == numSpecies, na.rm=TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityLight] <- sum(E$Status == 4 & E$SpeciesID == numSpecies, na.rm=TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityCompetition] <- sum(E$Status == 2 & E$SpeciesID == numSpecies, na.rm=TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityNatural] <- sum(E$Status == 5 & E$SpeciesID == numSpecies, na.rm=TRUE)
                # Climate mortality
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityHum] <- sum(E$Status == 6 & E$SpeciesID == numSpecies, na.rm=TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityTemp] <- sum(E$Status == 7 & E$SpeciesID == numSpecies, na.rm=TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityWind] <- sum(E$Status == 8 & E$SpeciesID == numSpecies, na.rm=TRUE)

                if (sum(E$Status == 1 & E$SpeciesID == numSpecies, na.rm=TRUE) > 0 && IntialNumberIndividuals[numSpecies] > 0) {
                    SummaryMatrixSpecies[rowIndex, ColSNumberPopulationGrowthRate] <- SummaryMatrixSpecies[rowIndex, ColSNumberIndividualsEnd] / SummaryMatrixSpecies[rowIndex, ColSNumberIndividualsBeginning]
                    SummaryMatrixSpecies[rowIndex, ColSNumberPopulationGrowthRateLog] <- log(SummaryMatrixSpecies[rowIndex, ColSNumberPopulationGrowthRate])
                    SummaryMatrixSpecies[rowIndex, ColSNumberBirthRate] <- NumberRecruitsPerSpecies[numSpecies] / IntialNumberIndividuals[numSpecies]
                    death_statuses <- c(2, 3, 4, 5) # 2: competition, 3: branch fall, 4: light, 5: natural mortality, 6: humidity, 7: temperature, 8: wind
                    SummaryMatrixSpecies[rowIndex, ColSNumberDeathRate] <- sum(E$Status %in% death_statuses & E$SpeciesID == numSpecies, na.rm = TRUE) / IntialNumberIndividuals[numSpecies]
                    SummaryMatrixSpecies[rowIndex, ColSAverageSize] <- mean(E$Mass[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSAverageAge] <- mean(E$Age[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMinLight] <- min(E$LightInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMaxLight] <- max(E$LightInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMeanLight] <- mean(E$LightInVoxel[E$SpeciesID == numSpecies])
                    # SummaryMatrixSpecies[rowIndex, ColSMinHeight] <- min(E$Z[E$SpeciesID == numSpecies])
                    # SummaryMatrixSpecies[rowIndex, ColSMaxHeight] <- max(E$Z[E$SpeciesID == numSpecies])
                    # SummaryMatrixSpecies[rowIndex, ColSMeanHeight] <- mean(E$Z[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMinHum] <- min(E$HumInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMaxHum] <- max(E$HumInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMeanHum] <- mean(E$HumInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMinTemp] <- min(E$TempInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMaxTemp] <- max(E$TempInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMeanTemp] <- mean(E$TempInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMinWind] <- min(E$WindInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMaxWind] <- max(E$WindInVoxel[E$SpeciesID == numSpecies])
                    SummaryMatrixSpecies[rowIndex, ColSMeanWind] <- mean(E$WindInVoxel[E$SpeciesID == numSpecies])
                } else {
                    # List of columns to be set to NaN
                    cols4summary <- c(
                      ColSNumberPopulationGrowthRate,
                      ColSNumberPopulationGrowthRateLog,
                      ColSNumberBirthRate,
                      ColSNumberDeathRate,
                      ColSAverageSize,
                      ColSAverageAge,
                      ColSMinLight,
                      ColSMaxLight,
                      ColSMeanLight,
                      ColSMinHum,
                      ColSMaxHum,
                      ColSMeanHum,
                      ColSMinTemp,
                      ColSMaxTemp,
                      ColSMeanTemp,
                      ColSMinWind,
                      ColSMaxWind,
                      ColSMeanWind
                    )
                    # Set all specified columns to NaN for the given row
                    SummaryMatrixSpecies[rowIndex, cols4summary] <- NaN
                    # SummaryMatrixSpecies[rowIndex, ColSMinHeight] <- NaN
                    # SummaryMatrixSpecies[rowIndex, ColSMaxHeight] <- NaN
                    # SummaryMatrixSpecies[rowIndex, ColSMeanHeight] <- NaN
                }

                SummaryMatrixSpeciesSave[rowIndex, 1] <- InitialTimeStep + t - 1
                SummaryMatrixSpeciesSave[rowIndex, int_seq(2, TotalColsSpeciesMatrix + 1)] <- SummaryMatrixSpecies[rowIndex, ]
            }

            ###############################################################################
            # Store information in SummaryMatrixCommunity
            SummaryMatrixCommunity$timeStep[t] <- InitialTimeStep + t - 1  # TimeStep
            SummaryMatrixCommunity$NumberSpeciesBeginning[t] <- InitialNumberSpecies  # NumberOfSpecies at beginning
            SummaryMatrixCommunity$NumberSpeciesEnd[t] <- length(unique(E$SpeciesID[E$Status == 1]))  # NumberOfSpecies at end
            SummaryMatrixCommunity$NumberIndividualsBeginning[t] <- IntialNumberIndividualsTotal  # NumberIndividuals at beginning
            SummaryMatrixCommunity$NumberIndividualsEnd[t] <- length(which(E$Status == 1))  # NumberIndividuals at end
            SummaryMatrixCommunity$Recruits[t] <- NumberRecruits  # Recruits
            SummaryMatrixCommunity$MortalityBranchFall[t] <- MortalityBranchFall  # MortalityBranchFall
            SummaryMatrixCommunity$MortalityLight[t] <- MortalityLight  # MortalityLight
            SummaryMatrixCommunity$MortalityCompetition[t] <- MortalityCompetition  # MortalityCompetition
            SummaryMatrixCommunity$MortalityNatural[t] <- MortalityNatural  # MortalityNatural
            SummaryMatrixCommunity$MortalityHum[t] <- MortalityHum  # Mortality due to humidity
            SummaryMatrixCommunity$MortalityTemp[t] <- MortalityTemp  # Mortality due to temperature
            SummaryMatrixCommunity$MortalityWind[t] <- MortalityWind  # Mortality due to wind
            SummaryMatrixCommunity$BranchSurfaceIndex[t] <- sum(Microhabitat[, , , Inds["TotalSurfaceAreaOpt"]]) / (dimPlot[1] * dimPlot[2])  # BranchSurfaceIndex
            SummaryMatrixCommunity$EpiphyteFilling[t] <- (sum(E$Mass^(2/3)) / SurfaceBiomassScaling) / sum(Microhabitat[, , , Inds["TotalSurfaceAreaOpt"]])  # EpiphyteFilling
            ###############################################################################

            # Command window information
            information <- "--------------------------------------------"
            information <- paste(information, paste0("Species Pool: ", numPool), sep="\n")
            information <- paste(information, paste0("Replicate: ", r), sep="\n")
            information <- paste(information, paste0("Time step: ", InitialTimeStep + t - 1), sep="\n")
            information <- paste(information, paste0("Number of individuals: ", SummaryMatrixCommunity$NumberIndividualsEnd[t]), sep="\n")
            information <- paste(information, paste0("Number of species: ", SummaryMatrixCommunity$NumberSpeciesEnd[t]), sep="\n")
            information <- paste(information, paste0("Number of recruits: ", NumberRecruits), sep="\n")
            information <- paste(information, paste0("MortalityBranchFall: ", MortalityBranchFall), sep="\n")
            information <- paste(information, paste0("MortalityLight: ", MortalityLight), sep="\n")
            information <- paste(information, paste0("MortalityCompetition: ", MortalityCompetition), sep="\n")
            information <- paste(information, paste0("MortalityNatural: ", MortalityNatural), sep="\n")
            information <- paste(information, paste0("MortalityHumidity: ", MortalityHum), sep="\n")
            information <- paste(information, paste0("MortalityTemperature: ", MortalityTemp), sep="\n")
            information <- paste(information, paste0("MortalityWind: ", MortalityWind), sep="\n")
            information <- paste(information, paste0("Time: ", format(Sys.time(), "%H:%M:%OS3")), sep="\n")
            writeLines(information)
            ###############################################################################

            # Saving
            # Save Epiphyte matrix for every time step
            ColumsToSave <- c("SpeciesID", "IndividualID", "Status", "Mass", "Age", "X", "Y", "Z", "TotalSurfaceInVoxel", "SurfaceLossInVoxel", "LightInVoxel", "HumInVoxel", "TempInVoxel", "WindInVoxel")
            write.csv(E[, ColumsToSave], file.path(DirectoryModelResultsRun, paste0("IndividualMatrixTimeStep", InitialTimeStep + t - 1, ".csv")), row.names=FALSE)

            # Create dataframe from matrix (including headers)
            SummaryMatrixSpeciesSave_df <- as.data.frame(SummaryMatrixSpeciesSave)
            names(SummaryMatrixSpeciesSave_df) <- SummaryMatrixSpeciesHeaders

            # Save SummaryMatrixSpecies for every time step
            write.csv(SummaryMatrixSpeciesSave_df, file.path(DirectoryModelResultsRun, "SpeciesSummary.csv"), row.names=FALSE)

            # Save SummaryMatrixCommunity for every time step (overwrite old one)
            write.csv(SummaryMatrixCommunity, file.path(DirectoryModelResultsRun, "CommunitySummary.csv"), row.names=FALSE)
            ###############################################################################

            # Remove dead individuals from Epimatrix
            E <- E[E$Status <= 1, ]  # Remove rows where Status > 1
        }
        return(NULL)
    }
}


main()
