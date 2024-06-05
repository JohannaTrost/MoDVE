# Create microhabitat matrices

# ------------------- Parameters ----------------------- #
# Parameters that need to be specified/checked before running this script

# This parameter determines which type of microhatiat matrices are generated:
# 1: real GroIMP forest with dynamics
# 2: static GroIMP forest (only forest at timeStepStart is used)
MicrohabitatType <- 1

# Parameters of light model
kL <- 0.6  # light extinction coefficient
DistVoxToConsider <- 8  # How many ring around focal voxel to consider in light model (5 voxels in x and y direction)

# Choose the forest parameters that shall be calculated and stored in the microhabitat matrix
# (this list can be extended for possible new applications of the epiphyte model.
# 1: use this variable
# 0: do not use it
TotalSurfaceAreaOpt <- 1
SurfaceAreaLossOpt <- 1
LightConditionsOpt <- 1
AverageWeightedAngles <- 0

# Parameters that need to be specified if MicrohabitatType=1 or MicrohabitatType=2
# Directory of GroIMP files (this directory is stored in the Microhabitat folder so that the
# connection to the input GroIMP files is always clear)
DirectoryGroIMP <- "/PATH/TO/INPUT"
DirectorySaveMain <- "/PATH/TO/OUTPUT"
DirectorySaveFolder <- "MicrohabitatMat"

# Name under which the microhabitat matrices are saved (The name of the folder under
# which the microhabitat matrices are saved is standarized and only the name
# of the forest is required here, and only if MicrohabitatType=1 or MicrohabitatType=2)
NameForest <- "ForestModel_Best"
ReplicateForest <- 0

# start and end timestep
timeStepStart <- 1
timeStepEnd <- 40

# Parameters that need to be specified if MicrohabitatType=3
# The following parameters are only needed if MicrohabitatType=3
# Dimensions of the theoretical forest
ForestHeight <- 40
dimXTheoretical <- 50
dimyTheoretical <- 50
BAI <- 3  # branch area index for the static, theoretical forest
LAI <- 6  # leaf are index

# Additional parameters
# Names of essential GroIMP files
shootFile <- paste("shoots_replicate_", ReplicateForest, "_time_step_", sep="")
trunkFile <- paste("trees_replicate_", ReplicateForest, "_time_step_", sep="")
voxelFile <- paste("voxel_replicate_", ReplicateForest, "_time_step_", sep="")

C <- c(0, 0, 1)  # Vector orthogonal to plane of X and Y
TotalVoxels <- (DistVoxToConsider * 2 + 1) ^ 2  # Total number of adjacent voxels considered
MatrixDimension <- sum(c(TotalSurfaceAreaOpt, SurfaceAreaLossOpt, LightConditionsOpt, AverageWeightedAngles))

# Other parameters created during porting to R
model_dir_name <- "Model"
results_dir_name <- "Results"
forest_global_param_name <- "Forest_param_global.txt"
forest_pass_param_name <- paste("Forest_param_pass", ReplicateForest, ".txt", sep="")


# The following parameters are generated automatically
# Load dimensions of forest patch from global forest file
if (MicrohabitatType == 1 || MicrohabitatType == 2) {
    path_to_forest_global <- file.path(DirectoryGroIMP, model_dir_name, forest_global_param_name)
    GlobalForest <- read.table(path_to_forest_global, sep="\t", row.names=1)

    MaxX <- GlobalForest["MaxX", 1]
    MaxY <- GlobalForest["MaxY", 1]
    MaxZ <- GlobalForest["MaxZ", 1]
    dimPlot <- c(MaxX, MaxY, MaxZ)
    corridor <- GlobalForest["WidthCorridor", 1]
    dimX <- MaxX + 2 * corridor  # MaxX+2*Corridor
    dimY <- MaxY + 2 * corridor  # MaxY+2*Corridor
    dimZ <- MaxZ  # MaxZ
}


# Create folder to save the microhabitat matrices
# The names of the folders are standadized:
# MicrohabitatType=1: 'Microhabitat_NameOfForestModel_SpatialExtent_timeSteps'
# MicrohabitatType=2: 'Microhabitat_NameOfForestModel_SpatialExtent_timeStep'
# MicrohabitatType=3: 'Microhabitat_BAI_LAI_kL'
if (MicrohabitatType == 1) {
    NameMicrohabitatMatrix <- paste("Microhabitat_", NameForest, "_", dimPlot[1], "x", dimPlot[2], "x", dimPlot[3], "_Rep", ReplicateForest, sep="")
    forest_type_str <- "DynamicForests"
} else if (MicrohabitatType == 2) {
    NameMicrohabitatMatrix <- paste("Microhabitat_", NameForest, "_", dimPlot[1], "x", dimPlot[2], "x", dimPlot[3], "_Rep", ReplicateForest, sep="")
    forest_type_str <- "StaticForests"
} else if (MicrohabitatType==3) {
    NameMicrohabitatMatrix <- paste("Microhabitat_BAI", BAI, "_LAI", LAI, "_kL", kL, sep="")
    forest_type_str <- "UniformForests"
}

DirectoryMatrices <- file.path(DirectorySaveMain, forest_type_str, DirectorySaveFolder, NameMicrohabitatMatrix)
dir.create(DirectoryMatrices, recursive=TRUE)

# Copy global and pass forest file to microhabitat folder
if (MicrohabitatType == 1 || MicrohabitatType == 2) {
    file.copy(file.path(DirectoryGroIMP, model_dir_name, forest_global_param_name), DirectoryMatrices, overwrite=FALSE)
    file.copy(file.path(DirectoryGroIMP, model_dir_name, forest_pass_param_name), DirectoryMatrices, overwrite=FALSE)
}

