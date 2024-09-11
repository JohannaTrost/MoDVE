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
timeSteps <- 40  # Model for timeSteps beginning at the time step given by the initial distribution

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

            # Erzeugen der Distanzmatrix mit allen Distanzen zum
            centralPoint <- c(floor(dimX/2) + 1, floor(dimY/2) + 1, floor(dimZ/2) + 1)
            DistanceMatrix <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))

            for (i in 1:dimX) {
                for (j in 1:dimY) {
                    for (k in 1:dimZ) {
                        # DistanceMatrix[i, j, k] <- pdist(rbind(c(i, j, k), c(centralPoint[1], centralPoint[2], centralPoint[3])))
                    }
                }
            }

        }
    }
}
