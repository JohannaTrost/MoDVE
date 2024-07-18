# Create the initial epiphyte distrubution depending on the epiphyte traits and the initial microhabitat matrix

###############################################################################
# Parameters that need to be specified/checked before running this script

# Name of epiphyte model
# FolderEpiphyteModel='EM_20160213'; %I should think about naming
SingleSpeciesModel <- 0  # 1: Single species model, 0: Community model

# Name of species pool
FolderSpeciesPools <- "SP_Random_IntAgeMat_2_IntRec_70_TraitCorrOn"

# Directory where model is save and directory where microhabitat matrices
# are stored
DirectoryModelMain <- "path/to/output"
DirectoryMicrohabitatMain <- "path/to/microhabitat"
DirectorySpeciesPoolsMain <- "path/to/species"

# Name of microhabitat matrix
FolderMicrohabitat <- "ForestModel_Best_30x30"
Replicate <- 0
MicrohabitatType <- 1  # Define which type of forest the microhabitat belongs to. 1: dynamic forest, 2: static forest, 3: uniform forest

# Choose species pools to use and number of replicates per species pool
numSpeciesPools <- c(99, 100)  # Start and end number of  species pools
replicatePerSpeciesPool <- 1  # Number of replicates per species pool
TimeStep <- 1  # Time step for which the Initial distribution is generated

# The suitable voxel can either be the voxel
# with the highest available surface area (MethodVoxel=1)
# or a random voxel (MethodVoxel=0)
MethodVoxel <- 0

# Define how many individuals per species are used, and how many of them are initially mature
# This variable defines if the NumberSpecies are total numbers irrespective
# of the model area (ScalingPerHa=0), or if the NumberSpecies or given per
# hectar and are scaled to the model area (ScalingPerHa=1)
ScalingPerHa <- 0
IndividualsPerSpecies <- 100
PercentageMaturePerSpecies <- 50

# This parameter set the scaling between the
SurfaceBiomassScaling <- 10000 * 0.01  # cm^2 per m^2
Imax <- 900
###############################################################################
