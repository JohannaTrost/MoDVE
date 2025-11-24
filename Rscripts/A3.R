# Create the initial epiphyte distrubution depending on the epiphyte traits and the initial microhabitat matrix
source("R/utils.R")

# Parse input configuration file
config <- parse_config("tests/config_a3.toml")
list2env(config, envir = globalenv())

###############################################################################
# Parameters that need to be specified/checked before running this script
{
  # RNG seed
  seed <- config$seed
  set.seed(seed, kind="Mersenne-Twister")  # integer for fixed seed or NULL for random

  SingleSpeciesModel <- config$SingleSpeciesModel  # 1: Single species model, 0: Community model

  # Directory where model is save and directory where microhabitat matrices
  # are stored
  DirectoryModelMain <- config$DirectoryModelMain
  DirectoryMicrohabitatMain <- config$DirectoryMicrohabitat
  DirectorySpeciesPoolsMain <- config$DirectorySpeciesPools

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
}

###############################################################################
# Generate directories to save the model
dir.create(DirectoryModelMain, recursive=TRUE)

###############################################################################
# Load parameters saved along with the microhabitat and species pool files

# Load plot dimensions if an artifical theoretical forest is used
dimPlot <- readRDS(file.path(DirectoryMicrohabitatMain, "dimPlot.rds"))

# Calculate individuals per species if normalization per hectare (ScalingPerHa=1) is chosen
if (ScalingPerHa == 1) {
    IndividualsPerSpecies <- IndividualsPerSpecies *
      ((dimPlot[1] * dimPlot[2]) / 10000)
}

# Get number of species from species pool file
species_filename <- paste("SpeciesPool", numSpeciesPools[1], ".csv", sep="")
Input_file <- file.path(DirectorySpeciesPoolsMain, species_filename)
SpeciesPool <- read.csv(Input_file)
NumberSpecies <- length(SpeciesPool$SpeciesID)

ColumnHeaders <- c("X", "Y", "Z", "Mass", "Status", "IndividualID", "SurfaceAreaOccupied", "Age")

# Get numbers of columns used in this script
ColMinLight <- match("MinLight", SpeciesPool)
ColMaxLight <- match("MaxLight", SpciesPool)
col_x <- match("X", ColumnHeaders)
col_y <- match("Y", ColumnHeaders)
col_z <- match("Z", ColumnHeaders)
col_mass <- match("Mass", ColumnHeaders)
col_status <- match("Status", ColumnHeaders)
col_id <- match("IndividualID", ColumnHeaders)
col_sa <- match("SurfaceAreaOccupied", ColumnHeaders)
col_age <- match("Age", ColumnHeaders)

# Get number of total individuals for each replicate
TotalIndividuals <- NumberSpecies * IndividualsPerSpecies
NumberMaturesPerSpecies <- round(IndividualsPerSpecies *
                                   PercentageMaturePerSpecies / 100)

# Load initial microhabitat matrix
microhabitat_filename <- paste("MicrohabitatMatrix", TimeStep, ".rds", sep="")
FileInitalMatrix <- file.path(DirectoryMicrohabitatMain, microhabitat_filename)
Microhabitat <- readRDS(FileInitalMatrix)

# Set real light values (in microhabitat, relative light values are saved)
Microhabitat[, , , 3] <- Microhabitat[, , , 3] * Imax

###############################################################################

for (numPool in numSpeciesPools[1]:numSpeciesPools[2]) {

  print(paste("Number species pool:", numPool))

  # Load the species pool
  species_filename <- paste("SpeciesPool", numPool, ".csv", sep = "")
  Input_file <- file.path(DirectorySpeciesPoolsMain, species_filename)
  SpeciesPool <- read.csv(Input_file)

  SizeSpeciesPool <- dim(SpeciesPool)

  for (numReplicates in seq_len(replicatePerSpeciesPool)) {

    print(paste("Number replicate:", numReplicates))

    # Initialize epiphyte matrix
    init_ind_mat <- array(
      rep(0, TotalIndividuals * 8),
      dim = c(TotalIndividuals, 8)
    )

    # Initialize Available surface area
    AvailableSurfaceArea <- Microhabitat[, , , 1]  # Matrix to trace the still available surface area per voxel

    # Fill init_ind_mat with species trait informations and the
    # initial size of each individual
    for (numSpecies in seq_len(NumberSpecies)) {
      for (numIndividual in seq_len(IndividualsPerSpecies)) {

        idx1 <- (numSpecies - 1) * IndividualsPerSpecies + numIndividual

        # Get size of individual
        if (numIndividual <= NumberMaturesPerSpecies) {
          SizeOfIndividual <- runif(1, SpeciesPool$MassAtMaturity[numSpecies],
                                    SpeciesPool$MaximumMass[numSpecies])  # Size of mature individuals
        } else {
          SizeOfIndividual <- runif(1, min = 0, max = SpeciesPool$MassAtMaturity[numSpecies])  # Size of juvenile individuals
        }

        # Store initial size of individual
        init_ind_mat[idx1, col_mass] <- SizeOfIndividual

        # store initial age of individual (age when it would have grown under optimal conditions)
        init_ind_mat[idx1, col_age] <- round(
          -log(1 - (SizeOfIndividual / SpeciesPool$MaximumMass[numSpecies])) /
          SpeciesPool$GrowthRate[numSpecies]
          )
      }
    }

    # Store individual ID for each individual
    init_ind_mat[seq_len(TotalIndividuals), col_id] <- seq_len(TotalIndividuals)

    # Calculate the surface area needed to support an individual of
    # this size =SurfaceAreaNeededInVoxel
    init_ind_mat[, col_sa] <- (init_ind_mat[, col_mass]^(2 / 3)) / SurfaceBiomassScaling

    if (SingleSpeciesModel) {
      # loop over all species
      for (numSpecies in seq_len(NumberSpecies)) {

        # Get subset of individuals for each species
        init_ind_matSub <- init_ind_mat[init_ind_mat[, 1] == numSpecies, ]

        # loop through all individuals, beginning with the largest (competition)
        # sort init_ind_matSub by size
        sort_ids <- order(init_ind_matSub[, col_mass], decreasing = TRUE)
        init_ind_matSub <- init_ind_matSub[sort_ids, ]
        NumNoSurface <- 0

        # Calculate potential voxels for each species which fulfill
        # their niche requirements (to save time they are precomputed here)
        comp1 <- Microhabitat[, , , 1] > 0
        comp2 <- Microhabitat[, , , 3] >= SpeciesPool$MinLight[numSpecies]
        comp3 <- Microhabitat[, , , 3] <= SpeciesPool$MaxLight[numSpecies]
        eligibleVoxels <- arrayInd(which(comp1 & comp2 & comp3), dim(Microhabitat))
        xcoords <- eligibleVoxels[, 1]
        ycoords <- eligibleVoxels[, 2]
        zcoords <- eligibleVoxels[, 3]

        for (ind in seq_len(nrow(init_ind_matSub))) { # loop over individuals

          # Create a random permutation
          randVoxels <- sample(seq_len(length(xcoords)))

          for (v in randVoxels) {

            x <- xcoords[v]
            y <- ycoords[v]
            z <- zcoords[v]

            avail_surf_area <- AvailableSurfaceArea[x, y, z]
            reqd_surf_area <- init_ind_matSub[ind, col_sa]

            if (avail_surf_area > reqd_surf_area) {

              # Allocate individual to voxel
              init_ind_matSub[ind, col_x] <- x
              init_ind_matSub[ind, col_y] <- y
              init_ind_matSub[ind, col_z] <- z

              # Set status of individual: status=1 => alive
              init_ind_matSub[ind, col_status] <- 1

              # Allocated individual reduces available surface accordingly
              AvailableSurfaceArea[x, y, z] <- AvailableSurfaceArea[x, y, z] -
                reqd_surf_area

              break # no need to consider other voxels
            }
            # else, move on to the next potential voxel

          } # random voxel loop
        } # individual loop

        # Set the coordinates to 1 for all individuals that did not find any
        # suitable habitat
        init_ind_matSub[init_ind_matSub[, col_x] == 0, col_status] <- 2
        init_ind_matSub[init_ind_matSub[, col_x] == 0, col_x] <- 1
        init_ind_matSub[init_ind_matSub[, col_y] == 0, col_y] <- 1
        init_ind_matSub[init_ind_matSub[, col_z] == 0, col_z] <- 1

        # Re-allocate subset to main matrix
        init_ind_mat[int_seq(
          from = ((numSpecies - 1) * IndividualsPerSpecies) + 1,
          to = numSpecies * IndividualsPerSpecies,
          by = 1
        ), ] <- init_ind_matSub

      }

    } else { # community model

      # loop randomly through all individuals and select suitable
      # voxel for each.
      RandNumInd <- sample(seq_len(TotalIndividuals), TotalIndividuals, replace=FALSE)

      for (i in seq_len(TotalIndividuals)) {
        print(paste("Individual number", i))

        NumIndRand <- RandNumInd[i]

        # Find all suitable voxels for this individual
        MinLightInd <- init_ind_mat[NumIndRand, ColMinLight]
        MaxLightInd <- init_ind_mat[NumIndRand, ColMaxLight]
        AreaNeededInd <- init_ind_mat[NumIndRand, col_sa]

        # 1. Get the positions of all voxels fullfilling the
        # requirements of the individual (light+area)
        tmp1 <- AvailableSurfaceArea[, , ] > AreaNeededInd
        tmp2 <- Microhabitat[, , , 3] >= MinLightInd
        tmp3 <- Microhabitat[, , , 3] <= MaxLightInd
        SuitableVoxels <- which(tmp1 & tmp2 & tmp3)

        # Choose one of the suitable voxels based on the specified
        # Method (if suitable voxels are available)
        if (length(SuitableVoxels) > 0) {

          if (MethodVoxel == 1) { # voxel with the most available surface area
            whichVoxel <- which(AvailableSurfaceArea[SuitableVoxels] ==
                max(AvailableSurfaceArea[SuitableVoxels]))[1]
          } else { # random voxel
            whichVoxel <- sample.int(length(SuitableVoxels), size=1)
          }
          ids <- arrayInd(SuitableVoxels[whichVoxel], dim(AvailableSurfaceArea))
          x <- ids[, 1]
          y <- ids[, 2]
          z <- ids[, 3]

          # Update available Surface Area
          AvailableSurfaceArea[x, y, z] <- AvailableSurfaceArea[x, y, z] - AreaNeededInd

          # Update Initial Epiphyte Matrix
          init_ind_mat[NumIndRand, col_x] <- x
          init_ind_mat[NumIndRand, col_y] <- y
          init_ind_mat[NumIndRand, col_z] <- z

          # Set status of individual: status=1 => alive
          init_ind_mat[NumIndRand, col_status] <- 1

        } else { # no suitable voxels

          # Set status of individual: status=2 =>> dead
          init_ind_mat[NumIndRand, col_status] <- 2

          # Set coordinates to 1 (might cause problems in later model if not)
          init_ind_mat[NumIndRand, col_x] <- 1
          init_ind_mat[NumIndRand, col_y] <- 1
          init_ind_mat[NumIndRand, col_z] <- 1
        }
      } # loop individuals

    } # ifelse singleSpecies

    # Copy trait data from SpeciesPool to InitalEpiphyteMatrix
    init_ind_mat[idx1, seq_len(SizeSpeciesPool[2])] <- as.numeric(SpeciesPool[numSpecies, ])

    # Create dataframe from matrix (including headers)
    init_ind_mat_df <- as.data.frame(init_ind_mat)
    names(init_ind_mat_df) <- ColumnHeaders

    # Save Initial Epiphyte Matrix
    FileName <- paste(
      "ID_SpeciesP_", numPool, "_Rep_", numReplicates, ".csv", sep = ""
      )
    write.csv(
      init_ind_mat_df, file.path(DirectoryModelMain, FileName),row.names = FALSE
    )
  } # replicate loop
} # species pool loop


