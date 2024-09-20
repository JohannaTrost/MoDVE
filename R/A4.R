###############################################################################
# START DELETEME

# MicrohabitatNumber = 1
# InititalDistNumber = 1

# FolderEpiphyteModel = "ForestModel_Best_30x30_Rep0"
# FolderInitialDistribution = "SP_Random_IA_2_IR_70_TimeS_1"
# FolderInitialDistributionsTemp = "SP_Random_IA_2_IR_70_TimeS_1"
# InitialTimeStepTemp = "1"  # SP_Random_IA_2_IR_70_TimeS_1 <-- that 1
# InitialTimeStep = 1

# FolderEpiphyteModelMain = "DynamicForests"
# FolderModelType = "CommunityModels"

# DirectoryEpiphyteModel = "INPUT/CommunityModels/DynamicForests/ForestModel_Best_30x30_Rep0"
# DirectoryIntitalDistribution = [DirectoryEpiphyteModel]/IniDist/SP_Random_IA_2_IR_70_TimeS_1"

# END DELETEME
###############################################################################

# Suppress warnings
# options(warn=-1)


# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community


###############################################################################
# Parameters that need to be specified/checked before running this script

# Folder of epiphyte models (these models are simulated in this order)
# The names of the models in the Folder "EpiphyteModels" are needed here
FolderEpiphyteModels <- c("ForestModel_Best_30x30_Rep0")  # unused?
MicrohabitatType <- 1  # Define which type of forest the microhabitat belongs to. 1: dynamic forest, 2: static forest, 3: uniform forest
SingleSpeciesModel <- 0  # 1: Single species model, 0: Community model

# Choose initial distributions (have to be located in 'FolderEpiphyteModel\IniDist\')
FolderInitialDistributions <- c("SP_Random_IA_2_IR_70_TimeS_1")  # unused?

DirectoryModelMain <- "path/to/A3/output"

# Model parameters
timeSteps <- 39  # Model for timeSteps beginning at the time step given by the initial distribution

# Density of individuals per ha at which to stop the simulationof the community and
# move to the next replicate (to prevent exploding communities)
StopCriterionHa <- 3000000  # Individuals per ha

# Choose species pools to use and number of replicates per species pool
numSpeciesPools <- c(99, 100)  # Start and end number of  species pools (if the species pools do not exist, they are automatically skipped)
replicatePerSpeciesPool <- 1  # Number of replicates per species pool  (if the replicates do not exist, they are automatically skipped)

SurfaceBiomassScaling <- 100  # cm^2 per m^2
Imax <- 900  # maximum light above canopy

# Competition Methods; defines which individuals are removed in voxels which
# are entirely filled. 1:size (small individuals are outcompetet by larger ones); 2:random competition
CompetitionMethod <- 1

# Mortality method (complete random or scaling with mass according to metabolic theory);
MortalityMethod <- 1  # 0: random mortality; 1: scaling with mass to the exponent -1/4
MortRateRandom <- 0.1
MortRateMass <- 0.1
MortRateMassScaling <- -0.25  # widely used scaling fator


# NEW
# Directories where microhabitat and species pool is saved
DirectoryMicrohabitat <- "path/to/A1/output"  # array???
DirectorySpeciesPools <- "path/to/A2/output"  # array???
DirectoryModelResults <- "path/to/output"  # array???

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


for (MicrohabitatNumber in 1:length(FolderEpiphyteModels)) {
    for (InititalDistNumber in 1:length(FolderInitialDistributions)) {

        # Create folder to save the model results
        dir.create(DirectoryModelResults, recursive=TRUE)

        # Load plot dimensions
        dimPlot <- readRDS(file.path(DirectoryMicrohabitat, "dimPlot.rds"))

        # Set StopCriterion for this simulation
        StopCriterion <- StopCriterionHa * dimPlot[1] * dimPlot[2] * 0.0001

        # Load TraitRanges (ranges used to create the species pool)
        FileTraitRanges <- file.path(DirectorySpeciesPools, "TraitRanges.csv")
        TraitRanges <- read.table(FileTraitRanges, sep=",", header=FALSE)

        SlopeRecruitment <- TraitRanges[1, 1]
        InterceptRecruitment <- TraitRanges[2, 1]


        TotalColsSpeciesMatrix <- 22

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
                SummaryMatrixSpeciesSave <- array(rep(0, timeSteps*NumberOfSpecies * TotalColsSpeciesMatrix), dim=c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix))

                # Initialize Matrix where speceies parameters are save
                SummaryMatrixSpecies <- array(rep(0, timeSteps*NumberOfSpecies * TotalColsSpeciesMatrix), dim=c(timeSteps*NumberOfSpecies, TotalColsSpeciesMatrix))



            }

        }
    }
}
