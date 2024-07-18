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

###############################################################################
# Folders and directories (these files should not change)

# The maximum path name length can be exeeded thus the folder names are abbreviated
if (MicrohabitatType == 1) {
    FolderEpiphyteModel <- "DynamicForests"
} else if (MicrohabitatType == 2) {
    FolderEpiphyteModel <- "StaticForests"
} else if (MicrohabitatType == 3) {
    FolderEpiphyteModel <- "UniformForests"
}

# Create main model folder under which the initial distribution is saved
if (SingleSpeciesModel == 1) {
    DirectoryEpiphyteModelMain1 <- file.path(DirectoryModelMain, "SingleSpeciesModels")
    DirectoryEpiphyteModelMain2 <- file.path(DirectoryEpiphyteModelMain1, FolderEpiphyteModel)
    dirname <- paste(FolderMicrohabitat, "_Rep", Replicate, sep="")
    DirectoryEpiphyteModel <- file.path(DirectoryEpiphyteModelMain2, dirname)
} else if (SingleSpeciesModel == 0) {
    DirectoryEpiphyteModelMain1 <- file.path(DirectoryModelMain, "CommunityModels")
    DirectoryEpiphyteModelMain2 <- file.path(DirectoryEpiphyteModelMain1, FolderEpiphyteModel)
    dirname <- paste(FolderMicrohabitat, "_Rep", Replicate, sep="")
    DirectoryEpiphyteModel <- file.path(DirectoryEpiphyteModelMain2, dirname)
}

# Generate abbreveation for species pool (in IniDist folder)
PosUnderscores <- unlist(gregexpr("_", FolderSpeciesPools))
substr1 <- substr(FolderSpeciesPools, PosUnderscores[1] + 1, PosUnderscores[2] - 1)
substr2 <- substr(FolderSpeciesPools, PosUnderscores[3] + 1, PosUnderscores[4] - 1)
substr3 <- substr(FolderSpeciesPools, PosUnderscores[5] + 1, PosUnderscores[6] - 1)
NameSpeciesPoolSave <- paste("SP_", substr1, "_IA_", substr2, "_IR_", substr3, "_TimeS_", TimeStep, sep="")

DirectoryIntitalDistributionMain <- file.path(DirectoryEpiphyteModel, "IniDist")
DirectoryIntitalDistribution <- file.path(DirectoryIntitalDistributionMain, NameSpeciesPoolSave)

# Generate directories to save the model
dir.create(DirectoryIntitalDistribution, recursive=TRUE)

dirname <- paste("Microhabitat_", FolderMicrohabitat, "_Rep", Replicate, sep="")
if (MicrohabitatType == 1) {
    DirectoryMicrohabitat <- file.path(DirectoryMicrohabitatMain, "DynamicForests", FolderMicrohabitat, dirname)
} else if (MicrohabitatType == 2) {
    DirectoryMicrohabitat <- file.path(DirectoryMicrohabitatMain, "StaticForests", FolderMicrohabitat, dirname)
} else if (MicrohabitatType == 3) {
    DirectoryMicrohabitat <- file.path(DirectoryMicrohabitatMain, "UniformForests", FolderMicrohabitat, dirname)
}

DirectorySpeciesPools <- file.path(DirectorySpeciesPoolsMain, FolderSpeciesPools)

###############################################################################
# Load parameters saved along with the microhabitat and species pool files

# Load plot dimensions if an artifical theoretical forest is used
dimPlot <- readRDS(file.path(DirectoryMicrohabitatMain, "dimPlot.rds"))

# Calculate individuals per species if normalization per hectare (ScalingPerHa=1) is chosen
if (ScalingPerHa == 1) {
    IndividualsPerSpecies <- IndividualsPerSpecies * ((dimPlot[1] * dimPlot[2]) / 10000)
}

# Get number of species from species pool file
species_filename <- paste("SpeciesPool", numSpeciesPools[1], ".csv", sep="")
Input_file <- file.path(DirectorySpeciesPoolsMain, species_filename)
SpeciesPool <- read.csv(Input_file)
