# Suppress warnings
# options(warn=-1)


# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("utils.R")

###############################################################################

compute_prob_matrix_norm <- function(centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool) {
    # Erzeugen der Distanzmatrix mit allen Distanzen zum
    DistanceMatrix <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))

    for (i in 1:dimX) {
        for (j in 1:dimY) {
            for (k in 1:dimZ) {
                x1 <- c(i, j, k)
                x2 <- c(centralPoint[1], centralPoint[2], centralPoint[3])
                DistanceMatrix[i, j, k] <- sqrt(sum((x1 - x2)^2))  # call to pdist() in the matlab script
            }
        }
    }

    # Erzeugen der Wahrschienlichkeitsmatrix anhand der Distanzmatrix und dem
    # artspezischien Wert bb aus der Epiphytenmatrix
    # Google translate: Generating the probability matrix based on the distance
    # matrix and the species value from the epiphyte matrix
    # negExp = @(distance,bb) exp(-distance.*bb); %Negative Exponential function
    ProbabilityMatrix <- array(rep(0, dimX * dimY * dimZ * NumberOfSpecies), dim=c(dimX, dimY, dimZ, NumberOfSpecies))
    ProbabilityMatrixNormalized <- array(rep(0, dimX * dimY * dimZ * NumberOfSpecies), dim=c(dimX, dimY, dimZ, NumberOfSpecies))

    for (i in 1:NumberOfSpecies) {
        exponentE <- SpeciesPool$DispersalKernel[i]
        dispersalAsymmetry <- SpeciesPool$DispersalKernelAsymmetry[i]

        ProbabilityMatrix[, , , i] <- exp(-DistanceMatrix * exponentE)  # call to negExp(DistanceMatrix(:,:,:),exponentE) in matlab

        # Apply dispersal asymmetry (probability to disperse downwards higher than upwards dispersal)
        # WARNING: The index i in the 4th dimension was missing from the Matlab script and so we were
        # getting "incorrect number of dimensions" in R. Added it here but need to ask if this is what
        # it was supposed to do.
        ProbabilityMatrix[, , centralPoint[3]:dimZ, i] <- ProbabilityMatrix[, , centralPoint[3]:dimZ, i] * ((1 - dispersalAsymmetry) / 0.5)
        ProbabilityMatrix[, , 1:(centralPoint[3] - 1), i] <- ProbabilityMatrix[, , 1:(centralPoint[3] - 1), i] * (dispersalAsymmetry / 0.5)

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


main <- function() {
    # Parse input configuration file
    config <- parse_config()

    ###############################################################################
    # Parameters that need to be specified/checked before running this script

    # RNG seed
    seed <- config$seed
    set.seed(seed, kind="Mersenne-Twister")

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

    # Mortality method (complete random or scaling with mass according to metabolic theory);
    MortalityMethod <- config$MortalityMethod  # 0: random mortality; 1: scaling with mass to the exponent -1/4
    MortRateRandom <- config$MortRateRandom
    MortRateMass <- config$MortRateMass
    MortRateMassScaling <- config$MortRateMassScaling  # widely used scaling fator

    InitialTimeStep <- config$InitialTimeStep  # Time step for which the Initial distribution is generated in A3

    ###############################################################################

    # Create folder to save the model results
    dir.create(DirectoryModelResults, recursive=TRUE)

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
    ColSMinHeight <- 20
    ColSMaxHeight <- 21
    ColSMeanHeight <- 22
    TotalColsSpeciesMatrix <- ColSMeanHeight

    # Headers of matrix
    SummaryMatrixSpeciesHeaders <- c("TimeStep", "SpeciesID", "NumberIndividualsBeginning", "NumberIndividualsEnd",
        "NumberMatureIndividuals", "NumberRecruits", "NumberRecruitsPotential", "NumberMortalityBranchFall",
        "NumberMortalityLight", "NumberMortalityCompetition", "NumberMortalityNatural", "PopulationGrowthRate",
        "PopulationGrowthRateLog", "BirthRate", "DeathRate", "AverageMass", "AverageAge", "MinLight", "MaxLight",
        "MeanLight", "MinHeight", "MaxHeight", "MeanHeight")

    if (length(SummaryMatrixSpeciesHeaders) != (TotalColsSpeciesMatrix + 1)) {
        stop("Headers of species matrix do not match with number of columns")
    }

    # Headers of matrix
    SummaryMatrixCommunityHeaders <- c("timeStep", "NumberSpeciesBeginning", "NumberSpeciesEnd",
        "NumberIndividualsBeginning", "NumberIndividualsEnd", "Recruits", "MortalityBranchFall",
        "MortalityLight", "MortalityCompetition", "MortalityNatural", "BranchSurfaceIndex", "EpiphyteFilling")



    ###############################################################################
    # Main loop for the community model

    for (numPool in numSpeciesPools[1]:numSpeciesPools[2]) {

        # Check if a initial distribution for the species pool exists. If not, move on to the next species pool
        FileNameInitalDistributionPool <- paste("ID_SpeciesP_", numPool, "_Rep_1", ".csv", sep="")
        if (!file.exists(file.path(DirectoryModelMain, FileNameInitalDistributionPool))) {
            break
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

        ProbabilityMatrixNormalized <- compute_prob_matrix_norm(centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool)

        # Main model loop for each replicate
        for (r in 1:replicatePerSpeciesPool) {
            # Check if a initial distribution for the species pool exists. If not, move on to the next species pool
            FileNameInitalDistribution <- file.path(DirectoryModelMain, paste("ID_SpeciesP_", numPool, "_Rep_", r, ".csv", sep=""))
            if (!file.exists(FileNameInitalDistribution)) {
                break
            }

            # Create Save-Directory for each each replicate/initialDistribution
            DirectoryModelResultsRun <- file.path(DirectoryModelResults, paste("ID_SpeciesP_", numPool, "_Rep_", r, sep=""))
            dir.create(DirectoryModelResultsRun, recursive=TRUE)

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

            for (t in seq_len(timeSteps)) {

                # Check if the stop criterion is met
                if (length(which(E$Status == 1)) > StopCriterion) {
                    break
                }

                ###############################################################################
                # 1. Dispersal

                # Load microhabitat matrix for specific timeStep if dynamic forest is simulated
                if (MicrohabitatType == 1) {
                    Microhabitat <- readRDS(file.path(DirectoryMicrohabitat, paste("MicrohabitatMatrix", InitialTimeStep + t - 1, ".rds", sep="")))
                    Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax  # In the microhabitat matrix, the realtive light extinction is stored: convert to light values in ?mol*m-2*s-1
                    d1 <- dim(Microhabitat)[1]
                    d2 <- dim(Microhabitat)[2]
                    d3 <- dim(Microhabitat)[3]
                    pot_habitat <- array(rep(0, d1 * d2 * d3), dim=c(d1, d2, d3))
                }

                # Store number of individuals at beginning of time step
                # NOTE: We can probably move this outside the timestep loop
                # We need to first check whether the columns from E that we use
                # change or not.
                IntialNumberIndividuals <- array(rep(0, NumberOfSpecies))
                for (g in 1:NumberOfSpecies) {
                    # Count indices where SpeciesID is g and Status is 1
                    IntialNumberIndividuals[g] <- length(which(E$SpeciesID == g & E$Status == 1))
                }
                IntialNumberIndividualsTotal <- length(which(E$Status == 1))
                InitialNumberSpecies <- length(unique(E$SpeciesID[E$Status == 1]))
                NumberRecruitsPerSpecies <- array(rep(0, NumberOfSpecies))

                # Calculate free surface area per voxel
                AvailableSurfaceArea <- Microhabitat[, , , 1]
                for (i in seq_len(nrow(E))) {
                    SurfaceAreaNeededInVoxel <- E$Mass[i]^(2/3) / SurfaceBiomassScaling
                    AvailableSurfaceArea[E$X[i], E$Y[i], E$Z[i]] <- max(0, AvailableSurfaceArea[E$X[i], E$Y[i], E$Z[i]] - SurfaceAreaNeededInVoxel)
                }

                # Check if there are species left (~isempty(E) in matlab)
                if (nrow(E) > 0) {
                    unique_species <- unique(E$SpeciesID)  # list with species IDs of all present species

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

                            # Store potential normalized number of recruits in SummaryMatrixSpecies
                            SummaryMatrixSpecies[((i-1) * timeSteps) + t, ColSNumberRecruitsPotential] <- sum(ProbabilityMatrixPerSpecies)  # potential recruitment / sum(sum(sum(ProbabilityMatrixPerSpecies))) in matlab

                            # Matix containing all voxel for which the light requirements are fulfilled
                            # We use the first row from MatureIndividulsPerSpecies. Since its elements have the same SpeciesID
                            # then the MinLight and MaxLight is the same for all rows.
                            pot_habitat <- ifelse((Microhabitat[, , , 3] >= MatureIndividulsPerSpecies$MinLight[1]) & (Microhabitat[, , , 3] <= MatureIndividulsPerSpecies$MaxLight[1]), 1, 0)

                            # Final probabiliy matrix for new recruits
                            probability_recruits <- ProbabilityMatrixPerSpecies * pot_habitat * AvailableSurfaceArea

                            # Calculate number of recuits based on final probability matrix
                            Recruits <- array(rpois(length(probability_recruits), probability_recruits), dim=dim(probability_recruits))  # poissrnd(probability_recruits) in matlab

                            # Add new recruits to epiphyte matrix
                            num_recruits <- sum(Recruits)  # sum(sum(sum(Recruits))) in matlab
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
                }

                NumberRecruits <- length(which(E$Status == 1)) - IntialNumberIndividualsTotal
                ###############################################################################

                # Unclear what this line in the Matlab script is supposed to do.
                # From what I understand, the first column in E ("SpeciesID") takes non-zero values
                # only, so I think that E(:,1)==0 will always be empty.
                # E(E(:,1)==0,:)=[]; %in rare case, some individuals with only zeros are creates, which is wrong. This is to prevent the script to stop.

                ###############################################################################
                # Growth
                for (i in seq_len(nrow(E))) {
                    # maybe it is faster if I do not use the if statement => speed testing
                    if (E$Status[i] == 1) {
                        tmp1 <- GrowthRate(E$MaximumMass[i], E$Mass[i], E$GrowthRate[i])
                        tmp2 <- Parabol(E$LightResponseA[i], E$LightResponseB[i], E$LightResponseC[i], Microhabitat[E$X[i], E$Y[i], E$Z[i], 3])
                        E$Mass[i] <- E$Mass[i] + max(0, tmp1 * tmp2)
                    }

                    # Add information about the voxel to the epiphyte matrix
                    E$SurfaceAreaOccupied[i] <- (E$Mass[i]^(2/3)) / SurfaceBiomassScaling
                    E$TotalSurfaceInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 1]  # Total surface in voxel
                    E$SurfaceLossInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 2]  # Percentage surface loss in this year
                    E$LightInVoxel[i] <- Microhabitat[E$X[i], E$Y[i], E$Z[i], 3]  # Light conditions in voxel
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
                        if (!is.nan(Microhabitat[E$X[i], E$Y[i], E$Z[i], 2]) && runif(1, min=0, max=1) < Microhabitat[E$X[i], E$Y[i], E$Z[i], 2]) {  # Mortality due to branch fall
                            E$Status[i] <- 3
                        } else if (Microhabitat[E$X[i], E$Y[i], E$Z[i], 3] < E$MinLight[i] | Microhabitat[E$X[i], E$Y[i], E$Z[i], 3] > E$MaxLight[i]) {  # Mortality due to changing light conditions
                            E$Status[i] <- 4
                        } else if (MortalityMethod == 0 && runif(1, min=0, max=1) < MortRateRandom) {  # Natural mortality rate
                            E$Status[i] <- 5
                        } else if (MortalityMethod == 1 && runif(1, min=0, max=1) < (MortRateMass * (E$Mass[i]^MortRateMassScaling))) {
                            E$Status[i] <- 5
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
                        E[is.element(E$IndividualID, EpisInVoxel[(CumulativeSumOfSurfaceSum + 1):nrow(EpisInVoxel), "IndividualID"]), "Status"] <- 2
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

                ###############################################################################

                # Store information in SummaryMatrixSpecies (summary over time for each species
                for (numSpecies in seq_len(NumberOfSpecies)) {
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSSpeciesID] <- numSpecies
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberIndividualsBeginning] <- IntialNumberIndividuals[numSpecies]
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberIndividualsEnd] <- sum(E$Status == 1 & E$SpeciesID == numSpecies, na.rm=TRUE)
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberMatureIndividuals] <- sum(E$Status == 1 & E$SpeciesID == numSpecies & E$Mass >= E$MassAtMaturity, na.rm=TRUE)
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberRecruits] <- NumberRecruitsPerSpecies[numSpecies]
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberMortalityBranchFall] <- sum(E$Status == 3 & E$SpeciesID == numSpecies, na.rm=TRUE)
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberMortalityLight] <- sum(E$Status == 4 & E$SpeciesID == numSpecies, na.rm=TRUE)
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberMortalityCompetition] <- sum(E$Status == 2 & E$SpeciesID == numSpecies, na.rm=TRUE)
                    SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberMortalityNatural] <- sum(E$Status == 5 & E$SpeciesID == numSpecies, na.rm=TRUE)

                    if (sum(E$Status == 1 & E$SpeciesID == numSpecies, na.rm=TRUE) > 0 && IntialNumberIndividuals[numSpecies] > 0) {
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberPopulationGrowthRate] <- SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberIndividualsEnd] / SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberIndividualsBeginning]
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberPopulationGrowthRateLog] <- log(SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberPopulationGrowthRate])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberBirthRate] <- NumberRecruitsPerSpecies[numSpecies] / IntialNumberIndividuals[numSpecies]
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberDeathRate] <- (sum(E$Status == 3 & E$SpeciesID == numSpecies, na.rm=TRUE) + sum(E$Status == 4 & E$SpeciesID == numSpecies, na.rm=TRUE) + sum(E$Status == 2 & E$SpeciesID == numSpecies, na.rm=TRUE) + sum(E$Status == 5 & E$SpeciesID == numSpecies, na.rm=TRUE)) / IntialNumberIndividuals[numSpecies]
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSAverageSize] <- mean(E$Mass[E$SpeciesID == numSpecies])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSAverageAge] <- mean(E$Age[E$SpeciesID == numSpecies])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMinLight] <- min(E$LightInVoxel[E$SpeciesID == numSpecies])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMaxLight] <- max(E$LightInVoxel[E$SpeciesID == numSpecies])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMeanLight] <- mean(E$LightInVoxel[E$SpeciesID == numSpecies])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMinHeight] <- min(E$Z[E$SpeciesID == numSpecies])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMaxHeight] <- max(E$Z[E$SpeciesID == numSpecies])
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMeanHeight] <- mean(E$Z[E$SpeciesID == numSpecies])
                    } else {
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberPopulationGrowthRate] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberPopulationGrowthRateLog] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberBirthRate] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSNumberDeathRate] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSAverageSize] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSAverageAge] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMinLight] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMaxLight] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMeanLight] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMinHeight] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMaxHeight] <- NaN
                        SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ColSMeanHeight] <- NaN
                    }

                    SummaryMatrixSpeciesSave[((numSpecies-1) * timeSteps) + t, 1] <- InitialTimeStep + t - 1
                    SummaryMatrixSpeciesSave[((numSpecies-1) * timeSteps) + t, 2:(TotalColsSpeciesMatrix + 1)] <- SummaryMatrixSpecies[((numSpecies-1) * timeSteps) + t, ]
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
                SummaryMatrixCommunity$BranchSurfaceIndex[t] <- sum(Microhabitat[, , , 1]) / (dimPlot[1] * dimPlot[2])  # BranchSurfaceIndex
                SummaryMatrixCommunity$EpiphyteFilling[t] <- (sum(E$Mass^(2/3)) / SurfaceBiomassScaling) / sum(Microhabitat[, , , 1])  # EpiphyteFilling
                ###############################################################################

                # Command window information
                print("--------------------------------------------")
                print(paste("Time step", InitialTimeStep + t - 1, sep=" "))
                print(paste("Number of individuals", SummaryMatrixCommunity$NumberIndividualsEnd[t], sep=" "))
                print(paste("Number of species", SummaryMatrixCommunity$NumberSpeciesEnd[t], sep=" "))
                print(paste("Number of recruits", NumberRecruits, sep=" "))
                print(paste("MortalityBranchFall", MortalityBranchFall, sep=" "))
                print(paste("MortalityLight", MortalityLight, sep=" "))
                print(paste("MortalityCompetition", MortalityCompetition, sep=" "))
                print(paste("MortalityNatural", MortalityNatural, sep=" "))
                ###############################################################################

                # Saving
                # Save Epiphyte matrix for every time step
                ColumsToSave <- c("SpeciesID", "IndividualID", "Status", "Mass", "Age", "X", "Y", "Z", "TotalSurfaceInVoxel", "SurfaceLossInVoxel", "LightInVoxel")
                write.csv(E[, ColumsToSave], file.path(DirectoryModelResultsRun, paste("IndividualMatrixTimeStep", InitialTimeStep + t - 1, ".csv", sep="")), row.names=FALSE)

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
        }
    }
}


main()
