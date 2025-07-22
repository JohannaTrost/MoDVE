options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

setwd("/home/jtrost_ext/MoDVE/R/")

# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("utils.R")

library("foreach")
library("doParallel")
library("doRNG")
library("readr")

###############################################################################

compute_prob_matrix_norm <- function(centralPoint, dimPlot, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool, WindSpeed = NULL) {
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

    # This is in case no wind dispersal was specified by the user
    if (!is.null(WindSpeed)) {
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
    } else {
        # No effect of wind on dispersal by setting wind speed to 0
        WindSpeedFull <- 0
    }

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

# Functions used in the model

# Bertalanffy Growth
GrowthRate <- function(MaxMass, Mass, K) {
    return(K * (MaxMass - Mass))
}

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
                      Inds,
                      EnvVarFlags) {
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
                # Compute potential habitat based on selected flags
                vars <- c("LightSuitable", "HumSuitable", "TempSuitable", "WindSuitable")
                pot_habitat <- Reduce("&", mget(vars[as.logical(EnvVarFlags)]))

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
        numCores <- detectCores()-1
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
    DirectoryEnvScores <- config$DirectoryEnvScores

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
    UseWindDispersal <- config$UseWindDispersal

    # Mortality method (complete random or scaling with mass according to metabolic theory);
    MortalityMethod <- config$MortalityMethod  # 0: random mortality; 1: scaling with mass to the exponent -1/4
    MortRateRandom <- config$MortRateRandom
    MortRateMass <- config$MortRateMass
    MortRateMassScaling <- config$MortRateMassScaling  # widely used scaling fator

    # Get environmental variables to account for here
    EnvVarFlags <- config$EnvVarFlags

    # Get microhabitat matrix variables - dynamic handling of selected variables
    MicrohabitatVariableFlags <- config$MicrohabitatVariableFlags
    MhVarNames <- c("TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt", "AverageWeightedAngles",
                    "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
    # Only keep active options
    ActiveOpts <- MhVarNames[as.logical(MicrohabitatVariableFlags)]
    # Assign indices
    Inds <- setNames(seq_along(ActiveOpts), ActiveOpts)

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
                       .export=c("compute_prob_matrix_norm", "int_seq", "dispersal", "GrowthRate", "Parabol", "SuitabilityScore")) %dorng% {
        numPool <- pairs$numPool[pair_idx]
        r <- pairs$r[pair_idx]

        # Check if a initial distribution for the species pool exists. If not, move on to the next species pool
        FileNameInitalDistribution <- file.path(DirectoryModelMain, paste("ID_SpeciesP_", numPool, "_Rep_", r, ".csv", sep=""))
        if (!file.exists(FileNameInitalDistribution)) {
            return(NULL)
        }

        # First step: create probability matrices for each species
        # Load species pool
        SpeciesPoolFileName <- paste("SpeciesPool", numPool, ".csv", sep="")
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

        if (!UseWindDispersal) {
            ProbabilityMatrixNormalized <- compute_prob_matrix_norm(
                  centralPoint, dimPlot, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool
            )
        }

        # Load initial epiphyte distribution
        E <- read.csv(FileNameInitalDistribution, sep=",", header=TRUE)  # E for epiphytes

        # Add column to E for additional information
        E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0

        MaxIndividualID <- nrow(E)  # to trace individual IDs

        # Initialize Matrix where community parameters are save
        SummaryMatrixCommunity <- data.frame(matrix(0.0, nrow=timeSteps, ncol=length(SummaryMatrixCommunityHeaders)))
        colnames(SummaryMatrixCommunity) <- SummaryMatrixCommunityHeaders

        # Load microhabitat matrix if a uniform or static forest is simulated (only needs to be loaded once an not envery timestep)
        if (MicrohabitatType == 2 || MicrohabitatType == 3) {
            Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, "MicrohabitatMatrix1.rds"))
            Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax  # In the microhabitat matrix, the realtive light extinction is stored: convert to light values in ?mol*m-2*s-1

            d1 <- dim(Microhabitat)[1]
            d2 <- dim(Microhabitat)[2]
            d3 <- dim(Microhabitat)[3]
            pot_habitat <- array(rep(0, d1 * d2 * d3), dim=c(d1, d2, d3))
        }

        # Initialize matrices where the aggregated information on species level are saved
        SummaryMatrixSpeciesSave <- array(rep(0, (timeSteps*NumberOfSpecies) * (TotalColsSpeciesMatrix + 1)), dim=c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix + 1))

        # Initialize Matrix where speceies parameters are save
        SummaryMatrixSpecies <- array(rep(0, (timeSteps*NumberOfSpecies) * TotalColsSpeciesMatrix), dim=c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix))

        for (t in 0:timeSteps) {
            timeStep <- InitialTimeStep + t

            # Check if the stop criterion is met
            if (length(which(E$Status == 1)) > StopCriterion) {
                break
            }

            # Load microhabitat matrix for specific timeStep if dynamic forest is simulated
            Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, paste("MicrohabitatMatrix", timeStep, ".rds", sep="")))
            Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax  # In the microhabitat matrix, the realtive light extinction is stored: convert to light values in ?mol*m-2*s-1
            d1 <- dim(Microhabitat)[1]
            d2 <- dim(Microhabitat)[2]
            d3 <- dim(Microhabitat)[3]
            pot_habitat <- array(rep(0, d1 * d2 * d3), dim=c(d1, d2, d3))

            # Load environmental scores for the current time step
            envSuitPath <- file.path(DirectoryEnvScores,
                                     paste0("ScaledSuitability_", numPool, "_TimeStep", timeStep, ".h5"))
            EnvSuitScors <- rhdf5::h5read(envSuitPath, "ScaledSuitabilityScores")

            ###############################################################################
            # 1. Dispersal
            # Create probability matrix for each species
            if (UseWindDispersal) {
                ProbabilityMatrixNormalized <- compute_prob_matrix_norm(
                  centralPoint, dimPlot, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool,
                  Microhabitat[, , , Inds["WindNicheOpt"]])
            }

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
                Inds,
                EnvVarFlags
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
            all_suits <- c()
            all_suits_prec <- c()
            for (i in seq_len(nrow(E))) {
                # maybe it is faster if I do not use the if statement => speed testing
                if (E$Status[i] == 1) {
                    tmp1 <- GrowthRate(E$MaximumMass[i], E$Mass[i], E$GrowthRate[i])
                    tmp2 <- EnvSuitScors[E$X[i], E$Y[i], E$Z[i], E$SpeciesID[i]]

                    if (is.na(tmp2) | is.nan(tmp2)) {
                        tmp2 <- 0  # If the suitability is NA, set it to 0
                        message("Suitability score is NA.")
                    }

                    all_suits_prec <- c(all_suits_prec, tmp2)
                    MassGained <- max(0, tmp1 * tmp2)
                    E$Mass[i] <- E$Mass[i] + MassGained
                }

                # Add information about the voxel to the epiphyte matrix
                E$SurfaceAreaOccupied[i] <- (E$Mass[i]^(2/3)) / SurfaceBiomassScaling
                E$TotalSurfaceInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 1]  # Total surface in voxel
                E$SurfaceLossInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 2]  # Percentage surface loss in this year
                E$LightInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["LightNicheOpt"]]  # Light conditions in voxel
                E$HumInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["HumNicheOpt"]]  # Humidity in voxel
                E$TempInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["TempNicheOpt"]]  # Temperature in voxel
                E$WindInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["WindNicheOpt"]]  # Wind in voxel
            }

            # Print average and std of suitability and response
            print(paste("Prec. Average suitability:", mean(all_suits_prec, na.rm=TRUE),
                        "Prec. median suitability:", median(all_suits_prec, na.rm=TRUE),
                        "Prec. Std suitability:", sd(all_suits_prec, na.rm=TRUE)))

            ###############################################################################

            ###############################################################################
            # Mortality
            vars <- c("LightSuitable", "HumSuitable", "TempSuitable", "WindSuitable")
            vars <- vars[as.logical(EnvVarFlags)]
            for (i in seq_len(nrow(E))) {
                if (E$Status[i] == 1) {

                    # The following comparison would fail without the is.nan check,
                    # because Microhabitat contains NaNs in some entries and
                    # in R a comparison with a NaN returns NA, not a boolean.
                    # Note: We call runif repeatedly intentionally. See Issue #16 on Github
                    if (!is.nan(Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["SurfaceAreaLossOpt"]]) &&
                      runif(1, min = 0, max = 1) < Microhabitat[E$X[i], E$Y[i], E$Z[i], Inds["SurfaceAreaLossOpt"]]) {  # Mortality due to branch fall
                        E$Status[i] <- 3
                    } else if ("LightSuitable" %in% vars && E$LightInVoxel[i] < E$MinLight[i] | E$LightInVoxel[i] > E$MaxLight[i]) {  # Mortality due to changing light conditions
                        E$Status[i] <- 4
                    } else if (MortalityMethod == 0 && runif(1, min = 0, max = 1) < MortRateRandom) {  # Natural mortality rate
                        E$Status[i] <- 5
                    } else if (MortalityMethod == 1 && runif(1, min = 0, max = 1) < (MortRateMass * (E$Mass[i]^MortRateMassScaling))) {
                        E$Status[i] <- 5
                    } else if ("HumSuitable" %in% vars && !is.na(E$HumInVoxel[i]) && (E$HumInVoxel[i] < E$MinHum[i] | E$HumInVoxel[i] > E$MaxHum[i])) {
                        E$Status[i] <- 6
                    } else if ("TempSuitable" %in% vars && !is.na(E$TempInVoxel[i]) && (E$TempInVoxel[i] < E$MinTemp[i] | E$TempInVoxel[i] > E$MaxTemp[i])) {
                        E$Status[i] <- 7
                    } else if ("WindSuitable" %in% vars && !is.na(E$WindInVoxel[i]) && E$WindInVoxel[i] > E$MaxWind[i]) {
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
            ind_tmp <- arrayInd(which(TotalSurfaceArePerVoxelOccupied > Microhabitat[, , , 1]), dim(TotalSurfaceArePerVoxelOccupied))
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
                CumulativeSumOfSurfaceSum <- length(which(CumulativeSumOfSurface <= Microhabitat[IndX[i], IndY[i], IndZ[i], 1]))

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
                SummaryMatrixSpecies[rowIndex, ColSNumberIndividualsEnd] <- sum(E$Status == 1 & E$SpeciesID == numSpecies, na.rm = TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMatureIndividuals] <- sum(E$Status == 1 &
                                                                                     E$SpeciesID == numSpecies &
                                                                                     E$Mass >= E$MassAtMaturity, na.rm = TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberRecruits] <- NumberRecruitsPerSpecies[numSpecies]
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityBranchFall] <- sum(E$Status == 3 & E$SpeciesID == numSpecies, na.rm = TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityLight] <- sum(E$Status == 4 & E$SpeciesID == numSpecies, na.rm = TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityCompetition] <- sum(E$Status == 2 & E$SpeciesID == numSpecies, na.rm = TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityNatural] <- sum(E$Status == 5 & E$SpeciesID == numSpecies, na.rm = TRUE)
                # Climate mortality
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityHum] <- sum(E$Status == 6 & E$SpeciesID == numSpecies, na.rm = TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityTemp] <- sum(E$Status == 7 & E$SpeciesID == numSpecies, na.rm = TRUE)
                SummaryMatrixSpecies[rowIndex, ColSNumberMortalityWind] <- sum(E$Status == 8 & E$SpeciesID == numSpecies, na.rm = TRUE)

                if (sum(E$Status == 1 & E$SpeciesID == numSpecies, na.rm = TRUE) > 0 && IntialNumberIndividuals[numSpecies] > 0) {
                    SummaryMatrixSpecies[rowIndex, ColSNumberPopulationGrowthRate] <- SummaryMatrixSpecies[rowIndex, ColSNumberIndividualsEnd] / SummaryMatrixSpecies[rowIndex, ColSNumberIndividualsBeginning]
                    SummaryMatrixSpecies[rowIndex, ColSNumberPopulationGrowthRateLog] <- log(SummaryMatrixSpecies[rowIndex, ColSNumberPopulationGrowthRate])
                    SummaryMatrixSpecies[rowIndex, ColSNumberBirthRate] <- NumberRecruitsPerSpecies[numSpecies] / IntialNumberIndividuals[numSpecies]
                    death_statuses <- c(2, 3, 4, 5, 6, 7, 8) # 2: competition, 3: branch fall, 4: light, 5: natural mortality, 6: humidity, 7: temperature, 8: wind
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

                SummaryMatrixSpeciesSave[rowIndex, 1] <- timeStep
                SummaryMatrixSpeciesSave[rowIndex, int_seq(2, TotalColsSpeciesMatrix + 1)] <- SummaryMatrixSpecies[rowIndex,]
            }

            ###############################################################################
            # Store information in SummaryMatrixCommunity
            SummaryMatrixCommunity$timeStep[t] <- timeStep  # TimeStep
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
            SummaryMatrixCommunity$BranchSurfaceIndex[t] <- sum(Microhabitat[, , , 1]) / (dimPlot[1] * dimPlot[2])  # BranchSurfaceIndex
            SummaryMatrixCommunity$EpiphyteFilling[t] <- (sum(E$Mass^(2/3)) / SurfaceBiomassScaling) / sum(Microhabitat[, , , 1])  # EpiphyteFilling
            ###############################################################################

            # Command window information
            information <- "--------------------------------------------"
            information <- paste(information, paste("Species Pool: ", numPool, sep=""), sep="\n")
            information <- paste(information, paste("Replicate: ", r, sep=""), sep="\n")
            information <- paste(information, paste("Time step: ", timeStep, sep=""), sep="\n")
            information <- paste(information, paste("Number of individuals: ", SummaryMatrixCommunity$NumberIndividualsEnd[t], sep=""), sep="\n")
            information <- paste(information, paste("Number of species: ", SummaryMatrixCommunity$NumberSpeciesEnd[t], sep=""), sep="\n")
            information <- paste(information, paste("Number of recruits: ", NumberRecruits, sep=""), sep="\n")
            information <- paste(information, paste("MortalityBranchFall: ", MortalityBranchFall, sep=""), sep="\n")
            information <- paste(information, paste("MortalityLight: ", MortalityLight, sep=""), sep="\n")
            information <- paste(information, paste("MortalityCompetition: ", MortalityCompetition, sep=""), sep="\n")
            information <- paste(information, paste("MortalityNatural: ", MortalityNatural, sep=""), sep="\n")
            information <- paste(information, paste0("MortalityHumidity: ", MortalityHum), sep = "\n")
            information <- paste(information, paste0("MortalityTemperature: ", MortalityTemp), sep = "\n")
            information <- paste(information, paste0("MortalityWind: ", MortalityWind), sep = "\n")
            information <- paste(information, paste("Time: ", format(Sys.time(), "%H:%M:%OS3"), sep=""), sep="\n")
            writeLines(information)
            ###############################################################################

            # Saving
            # Save Epiphyte matrix for every time step
            ColumsToSave <- c("SpeciesID", "IndividualID", "Status", "Mass", "Age", "X", "Y", "Z", "TotalSurfaceInVoxel", "SurfaceLossInVoxel", "LightInVoxel")
            write.csv(E[, ColumsToSave], file.path(DirectoryModelResultsRun, paste("IndividualMatrixTimeStep", timeStep, ".csv", sep="")), row.names=FALSE)

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