# Create the initial epiphyte distrubution depending on the epiphyte traits and the initial microhabitat matrix
source("utils.R")


# Growth function. Used here to approximate the age of the individuals
# MassFunctionOfAge=@(MaxMass,K,Age) (MaxMass*(1-exp(-K*(Age))));
AgeFunctionOfMass <- function(MaxMass, Mass, K) {
    return(-log(1 - (Mass / MaxMass)) / K)
}

# Parse input configuration file
config <- parse_config()

###############################################################################
# Parameters that need to be specified/checked before running this script

# RNG seed
seed <- config$seed
set.seed(seed, kind="Mersenne-Twister")  # integer for fixed seed or NULL for random

SingleSpeciesModel <- config$SingleSpeciesModel  # 1: Single species model, 0: Community model

# Directory where model is save and directory where microhabitat matrices
# are stored
DirectoryModelMain <- config$DirectoryModelMain
DirectoryMicrohabitatMain <- config$DirectoryMicrohabitatMain
DirectorySpeciesPoolsMain <- config$DirectorySpeciesPoolsMain

Replicate <- config$Replicate

# Choose species pools to use and number of replicates per species pool
numSpeciesPools <- config$numSpeciesPools  # Start and end number of  species pools
replicatePerSpeciesPool <- config$replicatePerSpeciesPool  # Number of replicates per species pool
TimeStep <- config$TimeStep  # Time step for which the Initial distribution is generated

# The suitable voxel can either be the voxel
# with the highest available surface area (MethodVoxel=1)
# or a random voxel (MethodVoxel=0)
MethodVoxel <- config$MethodVoxel

# Define how many individuals per species are used, and how many of them are initially mature
# This variable defines if the NumberSpecies are total numbers irrespective
# of the model area (ScalingPerHa=0), or if the NumberSpecies or given per
# hectar and are scaled to the model area (ScalingPerHa=1)
ScalingPerHa <- config$ScalingPerHa
IndividualsPerSpecies <- config$IndividualsPerSpecies
PercentageMaturePerSpecies <- config$PercentageMaturePerSpecies

# This parameter set the scaling between the
SurfaceBiomassScaling <- config$SurfaceBiomassScaling  # cm^2 per m^2
Imax <- config$Imax
###############################################################################

###############################################################################
# Generate directories to save the model
dir.create(DirectoryModelMain, recursive=TRUE)

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
NumberSpecies <- length(SpeciesPool$SpeciesID)

ColumnHeaders <- c(colnames(SpeciesPool), c("X", "Y", "Z", "Mass", "Status", "IndividualID", "SurfaceAreaOccupied", "Age"))

# Get numbers of columns used in this script
ColMinLight <- match("MinLight", ColumnHeaders)
ColMaxLight <- match("MaxLight", ColumnHeaders)

# Get number of total individuals for each replicate
TotalIndividuals <- NumberSpecies * IndividualsPerSpecies
NumberMaturesPerSpecies <- round(IndividualsPerSpecies * (PercentageMaturePerSpecies / 100))

# Load initial microhabitat matrix
microhabitat_filename <- paste("MicrohabitatMatrix", TimeStep, ".rds", sep="")
FileInitalMatrix <- file.path(DirectoryMicrohabitatMain, microhabitat_filename)
Microhabitat <- readRDS(FileInitalMatrix)

# Set real light values (in microhabitat, relative light values are saved)
Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax

###############################################################################
# Main loop for Single Species Model
if (SingleSpeciesModel == 1) {
# Start - Identical for both SingleSpeciesModel 0 and 1
    for (numPool in numSpeciesPools[1]:numSpeciesPools[2]) {
        print(paste("Number species pool:", numPool))

        # Laden des species pools
        species_filename <- paste("SpeciesPool", numPool, ".csv", sep="")
        Input_file <- file.path(DirectorySpeciesPoolsMain, species_filename)
        SpeciesPool <- read.csv(Input_file)

        SizeSpeciesPool <- dim(SpeciesPool)

        for (numReplicates in 1:replicatePerSpeciesPool) {
            print(paste("Number replicate:", numReplicates))

            # Initialize epiphyte matrix
            IntitalEpiphyteMatrix <- array(rep(0, TotalIndividuals * (SizeSpeciesPool[2] + 8)), dim=c(TotalIndividuals, SizeSpeciesPool[2] + 8))

            # Initialize Available surface area
            AvailableSurfaceArea <- Microhabitat[, , , 1]  # Matrix to trace the still available surface area per voxel

            # Fill InitialEpiphyteMatrix with species trait informations and the
            # initial size of each individual
            for (numSpecies in 1:NumberSpecies) {
                for (numIndividual in 1:IndividualsPerSpecies) {
                    idx1 <- ((numSpecies - 1) * IndividualsPerSpecies) + numIndividual

                    # Copy trait data from SpeciesPool to InitalEpiphyteMatrix
                    IntitalEpiphyteMatrix[idx1, 1:SizeSpeciesPool[2]] <- as.numeric(SpeciesPool[numSpecies, ])

                    # Get size of individual
                    if (numIndividual <= NumberMaturesPerSpecies) {
                        SizeOfIndividual <- runif(1, min=SpeciesPool$MassAtMaturity[numSpecies], max=SpeciesPool$MaximumMass[numSpecies])  # Size of mature individuals
                    } else {
                        SizeOfIndividual <- runif(1, min=0, max=SpeciesPool$MassAtMaturity[numSpecies])  # Size of juvenile individuals
                    }

                    # Store initial size of individual
                    IntitalEpiphyteMatrix[idx1, SizeSpeciesPool[2] + 4] <- SizeOfIndividual

                    # store initial age of individual (age when it would have grown under optimal conditions)
                    IntitalEpiphyteMatrix[idx1, SizeSpeciesPool[2] + 8] <- round(AgeFunctionOfMass(SpeciesPool$MaximumMass[numSpecies], SizeOfIndividual, SpeciesPool$GrowthRate[numSpecies]))
                }
            }

            # Store individual ID for each individual
            IntitalEpiphyteMatrix[1:TotalIndividuals, SizeSpeciesPool[2] + 6] <- 1:TotalIndividuals

            # Calculate the surface area needed to support an individual of
            # this size =SurfaceAreaNeededInVoxel
            IntitalEpiphyteMatrix[, SizeSpeciesPool[2] + 7] <- (IntitalEpiphyteMatrix[, SizeSpeciesPool[2] + 4]^(2 / 3)) / SurfaceBiomassScaling
# End - Identical for both SingleSpeciesModel 0 and 1

            # loop over all species
            for (numSpecies in 1:NumberSpecies) {
                # Get subset of indiduals for each species
                IntitalEpiphyteMatrixSub <- IntitalEpiphyteMatrix[IntitalEpiphyteMatrix[, 1] == numSpecies, ]

                # loop through all individuals, beginning with the largest (competition)
                # sort IntitalEpiphyteMatrixSub by size
                sort_ids <- order(IntitalEpiphyteMatrixSub[, SizeSpeciesPool[2] + 4], decreasing=TRUE)
                IntitalEpiphyteMatrixSub <- IntitalEpiphyteMatrixSub[sort_ids, ]
                NumNoSurface <- 0

                # Calculate potential voxels for for each species which fullfil
                # their niche requirments (to save time they are precomputed here)
                comp1 <- Microhabitat[, , , 1] > 0
                comp2 <- Microhabitat[, , , 3] >= SpeciesPool$MinLight[numSpecies]
                comp3 <- Microhabitat[, , , 3] <= SpeciesPool$MaxLight[numSpecies]
                ids <- arrayInd(which(comp1 & comp2 & comp3), dim(Microhabitat))
                x <- ids[, 1]
                y <- ids[, 2]
                z <- ids[, 3]
                w <- ids[, 4]

                for (i in 1:nrow(IntitalEpiphyteMatrixSub)) {
                    sequence <- seq_len(length(x))  # Works as expected even if x is empty

                    # Create a random permutation
                    randNumbers <- sample(sequence)

                    for (PotVoxels in sequence) {
                        avail_surf_area <- AvailableSurfaceArea[x[randNumbers[PotVoxels]], y[randNumbers[PotVoxels]], z[randNumbers[PotVoxels]]]
                        if (avail_surf_area > IntitalEpiphyteMatrixSub[i, SizeSpeciesPool[2] + 7]) {
                            IntitalEpiphyteMatrixSub[i, SizeSpeciesPool[2] + 1] <- x[randNumbers[PotVoxels]]
                            IntitalEpiphyteMatrixSub[i, SizeSpeciesPool[2] + 2] <- y[randNumbers[PotVoxels]]
                            IntitalEpiphyteMatrixSub[i, SizeSpeciesPool[2] + 3] <- z[randNumbers[PotVoxels]]

                            # Set status of individual: status=1 => alive
                            IntitalEpiphyteMatrixSub[i, SizeSpeciesPool[2] + 5] <- 1

                            AvailableSurfaceArea[x[randNumbers[PotVoxels]], y[randNumbers[PotVoxels]], z[randNumbers[PotVoxels]]] <- AvailableSurfaceArea[x[randNumbers[PotVoxels]], y[randNumbers[PotVoxels]], z[randNumbers[PotVoxels]]] - IntitalEpiphyteMatrixSub[i, SizeSpeciesPool[2] + 7]

                            break
                        }
                    }
                }
                # Set the coordinates to 1 for all individuals that did not find any
                # suitable habitat
                IntitalEpiphyteMatrixSub[IntitalEpiphyteMatrixSub[, SizeSpeciesPool[2] + 1] == 0, SizeSpeciesPool[2] + 5] <- 2
                IntitalEpiphyteMatrixSub[IntitalEpiphyteMatrixSub[, SizeSpeciesPool[2] + 1] == 0, SizeSpeciesPool[2] + 1] <- 1
                IntitalEpiphyteMatrixSub[IntitalEpiphyteMatrixSub[, SizeSpeciesPool[2] + 2] == 0, SizeSpeciesPool[2] + 2] <- 1
                IntitalEpiphyteMatrixSub[IntitalEpiphyteMatrixSub[, SizeSpeciesPool[2] + 3] == 0, SizeSpeciesPool[2] + 3] <- 1

                IntitalEpiphyteMatrix[(((numSpecies-1) * IndividualsPerSpecies) + 1):(numSpecies * IndividualsPerSpecies), ] <- IntitalEpiphyteMatrixSub
            }
            # Create dataframe from matrix (including headers)
            IntitalEpiphyteMatrix_df <- as.data.frame(IntitalEpiphyteMatrix)
            names(IntitalEpiphyteMatrix_df) <- ColumnHeaders

            # Save Inital Epiphyte Matrix
            FileName <- paste("ID_SpeciesP_", numPool, "_Rep_", numReplicates, ".csv", sep="")
            write.csv(IntitalEpiphyteMatrix_df, file.path(DirectoryModelMain, FileName), row.names=FALSE)
        }
    }
}

# Main loop for community model
if (SingleSpeciesModel == 0) {
# Start - Identical for both SingleSpeciesModel 0 and 1
    for (numPool in numSpeciesPools[1]:numSpeciesPools[2]) {
        print(paste("Number species pool:", numPool))

        # Laden des species pools
        species_filename <- paste("SpeciesPool", numPool, ".csv", sep="")
        Input_file <- file.path(DirectorySpeciesPoolsMain, species_filename)
        SpeciesPool <- read.csv(Input_file)

        SizeSpeciesPool <- dim(SpeciesPool)

        for (numReplicates in 1:replicatePerSpeciesPool) {
            print(paste("Number replicate:", numReplicates))

            # Initialize epiphyte matrix
            IntitalEpiphyteMatrix <- array(rep(0, TotalIndividuals * (SizeSpeciesPool[2] + 8)), dim=c(TotalIndividuals, SizeSpeciesPool[2] + 8))

            # Initialize Available surface area
            AvailableSurfaceArea <- Microhabitat[, , , 1]  # Matrix to trace the still available surface area per voxel

            # Fill InitialEpiphyteMatrix with species trait informations and the
            # initial size of each individual
            for (numSpecies in 1:NumberSpecies) {
                for (numIndividual in 1:IndividualsPerSpecies) {
                    idx1 <- ((numSpecies - 1) * IndividualsPerSpecies) + numIndividual

                    # Copy trait data from SpeciesPool to InitalEpiphyteMatrix
                    IntitalEpiphyteMatrix[idx1, 1:SizeSpeciesPool[2]] <- as.numeric(SpeciesPool[numSpecies, ])

                    # Get size of individual
                    if (numIndividual <= NumberMaturesPerSpecies) {
                        SizeOfIndividual <- runif(1, min=SpeciesPool$MassAtMaturity[numSpecies], max=SpeciesPool$MaximumMass[numSpecies])  # Size of mature individuals
                    } else {
                        SizeOfIndividual <- runif(1, min=0, max=SpeciesPool$MassAtMaturity[numSpecies])  # Size of juvenile individuals
                    }

                    # Store initial size of individual
                    IntitalEpiphyteMatrix[idx1, SizeSpeciesPool[2] + 4] <- SizeOfIndividual

                    # store initial age of individual (age when it would have grown under optimal conditions)
                    IntitalEpiphyteMatrix[idx1, SizeSpeciesPool[2] + 8] <- round(AgeFunctionOfMass(SpeciesPool$MaximumMass[numSpecies], SizeOfIndividual, SpeciesPool$GrowthRate[numSpecies]))
                }
            }

            # Store individual ID for each individual
            IntitalEpiphyteMatrix[1:TotalIndividuals, SizeSpeciesPool[2] + 6] <- 1:TotalIndividuals

            # Calculate the surface area needed to support an individual of
            # this size =SurfaceAreaNeededInVoxel
            IntitalEpiphyteMatrix[, SizeSpeciesPool[2] + 7] <- (IntitalEpiphyteMatrix[, SizeSpeciesPool[2] + 4]^(2 / 3)) / SurfaceBiomassScaling
# End - Identical for both SingleSpeciesModel 0 and 1

            # loop randomly through all individuals and select suitable
            # voxel for each.
            RandNumInd <- sample(seq_len(TotalIndividuals), TotalIndividuals, replace=FALSE)

            for (i in 1:TotalIndividuals) {
                print(paste("Individual number", i))

                NumIndRand <- RandNumInd[i]

                # Find all suitable voxels for this individual
                MinLightInd <- IntitalEpiphyteMatrix[NumIndRand, ColMinLight]
                MaxLightInd <- IntitalEpiphyteMatrix[NumIndRand, ColMaxLight]
                AreaNeededInd <- IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 7]

                # 1. Get the postions of all voxels fullfilling the
                # requirements of the individual (light+area)
                tmp1 <- AvailableSurfaceArea[, , ] > AreaNeededInd
                tmp2 <- Microhabitat[, , , 3] >= MinLightInd
                tmp3 <- Microhabitat[, , , 3] <= MaxLightInd
                SuitableVoxels <- which(tmp1 & tmp2 & tmp3)

                # Choose one of the suitable voxels based on the specified
                # Method (if suitable voxels are available)
                if (length(SuitableVoxels) > 0) {
                    if (MethodVoxel == 1) {
                        MaxVal <- which(AvailableSurfaceArea[SuitableVoxels] == max(AvailableSurfaceArea[SuitableVoxels]))
                        ids <- arrayInd(SuitableVoxels[MaxVal[1]], dim(AvailableSurfaceArea))
                        x <- ids[, 1]
                        y <- ids[, 2]
                        z <- ids[, 3]
                    } else if (MethodVoxel == 0) {
                        RandVal <- sample.int(length(SuitableVoxels), size=1)
                        ids <- arrayInd(SuitableVoxels[RandVal], dim(AvailableSurfaceArea))
                        x <- ids[, 1]
                        y <- ids[, 2]
                        z <- ids[, 3]
                    }

                    # Update available Surface Area
                    AvailableSurfaceArea[x, y, z] <- AvailableSurfaceArea[x, y, z] - AreaNeededInd

                    # Update Inital Epiphyte Matrix
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 1] <- x
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 2] <- y
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 3] <- z

                    # Set status of individual: status=1 => alive
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 5] <- 1
                } else {
                    # Set status of individual: status=2 =>> dead
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 5] <- 2

                    # Set coordinates to 1 (might cause problems in later model if not)
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 1] <- 1
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 2] <- 1
                    IntitalEpiphyteMatrix[NumIndRand, SizeSpeciesPool[2] + 3] <- 1
                }
            }
            # Create dataframe from matrix (including headers)
            IntitalEpiphyteMatrix_df <- as.data.frame(IntitalEpiphyteMatrix)
            names(IntitalEpiphyteMatrix_df) <- ColumnHeaders

            # Save Inital Epiphyte Matrix
            FileName <- paste("ID_SpeciesP_", numPool, "_Rep_", numReplicates, ".csv", sep="")
            write.csv(IntitalEpiphyteMatrix_df, file.path(DirectoryModelMain, FileName), row.names=FALSE)
        }
    }
}
