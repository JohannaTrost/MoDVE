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


# Generation of microhabitat matrix of static or dynamic forest (MicrohabitatType=1 or MicrohabitatType=2)
# Here, the choosen parameter (total surface, surface loss, light conditions,average angle) are calculated for each voxel in each timestep

if (MicrohabitatType == 1 || MicrohabitatType == 2) {

    # In a static forest, only the initial forest at time step timeStepStart is of interest
    if (MicrohabitatType == 2) {
        timeStepEnd <- timeStepStart + 1
    }

    for (i in timeStepStart:(timeStepEnd-1)) {
        print(paste("Time step", i))
        start_time <- Sys.time()

        # Load shoot and trunk files of actual and next timestep: Shoots at begin of year
        # and at the end of year/begin of next year
        src_dir <- file.path(DirectoryGroIMP, results_dir_name)

        shootFileNameOld <- paste(shootFile, i, ".txt", sep="")
        shootFileNameNew <- paste(shootFile, i+1, ".txt", sep="")

        trunkFileNameOld <- paste(trunkFile, i, ".txt", sep="")
        trunkFileNameNew <- paste(trunkFile, i+1, ".txt", sep="")

        if (i != timeStepStart) {
            ShootsBegin <- ShootsEnd
            ShootsEnd <- read.table(file.path(src_dir, shootFileNameNew), sep="\t", header=TRUE, skip=1)
            # ShootsEnd <- readr::read_tsv(file.path(src_dir, shootFileNameNew), skip=1, show_col_types=FALSE)

            TrunksBegin <- TrunksEnd
            TrunksEnd <- read.table(file.path(src_dir, trunkFileNameNew), sep="\t", header=TRUE, skip=8)
            # TrunksEnd <- readr::read_tsv(file.path(src_dir, trunkFileNameNew), skip=8, show_col_types=FALSE)
        } else {
            ShootsBegin <- read.table(file.path(src_dir, shootFileNameOld), sep="\t", header=TRUE, skip=1)
            ShootsEnd <- read.table(file.path(src_dir, shootFileNameNew), sep="\t", header=TRUE, skip=1)

            TrunksBegin <- read.table(file.path(src_dir, trunkFileNameOld), sep="\t", header=TRUE, skip=8)
            TrunksEnd <- read.table(file.path(src_dir, trunkFileNameNew), sep="\t", header=TRUE, skip=8)
        }

        # inititalize all matrices
        Mat_surface_per_cell <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))
        Mat_weighted_angle_per_cell <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))
        Mat_light_per_cell <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))
        Mat_surfaceloss_per_cell <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))
        Mat_leafArea_per_cell <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))

        # Get all branch segments that die during time step
        DeadSegments <- ShootsBegin$shootID[!is.element(ShootsBegin$shootID, ShootsEnd$shootID)]
        DeadSegments1 <- is.element(DeadSegments, ShootsBegin$shootID)
        locDeadSegments <- match(DeadSegments, ShootsBegin$shootID)  # nomatch=0?
        CounterDead <- 1
        TotalDead <- length(locDeadSegments)

        # Loop through all shoots and calculate total surface area, surface loss and weighted angles
        for (j in 1:nrow(ShootsBegin)) {

            # Get voxel the shoot is intersecting with
            # Here, we are using a simple approach because the shoots are
            # usually small only intersect with maximum two voxels
            UniqueX <- c(ceiling(ShootsBegin$xbegin[j]), ceiling(ShootsBegin$xend[j]))
            UniqueY <- c(ceiling(ShootsBegin$ybegin[j]), ceiling(ShootsBegin$yend[j]))
            UniqueZ <- c(ceiling(ShootsBegin$zbegin[j]), ceiling(ShootsBegin$zend[j]))

            numX <- sum((UniqueX[2] - UniqueX[1]) > 0, na.rm=TRUE) + 1
            numY <- sum((UniqueY[2] - UniqueY[1]) > 0, na.rm=TRUE) + 1
            numZ <- sum((UniqueZ[2] - UniqueZ[1]) > 0, na.rm=TRUE) + 1

            for (x in 1:numX) {
                xid <- UniqueX[x]

                for (y in 1:numY) {
                    yid <- UniqueY[y]

                    for (z in 1:numZ) {
                        zid <- UniqueZ[z]

                        # Calculate new total surface area per voxel
                        if (TotalSurfaceAreaOpt == 1) {
                            # Surface of single branch within voxel
                            voxel_branch_surface <- (ShootsBegin$length[j] / (numX*numY*numZ)) * ShootsBegin$diameter[j] * pi/2

                            Mat_surface_per_cell[xid, yid, zid] <- Mat_surface_per_cell[xid, yid, zid] + voxel_branch_surface
                        }

                        # If branch is lost during this time step, add it to lost surface
                        if (SurfaceAreaLossOpt == 1) {
                            if (j == locDeadSegments[CounterDead]) {
                                CounterDead <- min(TotalDead, CounterDead + 1)

                                lost_surface <- (ShootsBegin$length[j] / (numX*numY*numZ)) * ShootsBegin$diameter[j] * pi/2
                                Mat_surfaceloss_per_cell[xid, yid, zid] <- Mat_surfaceloss_per_cell[xid, yid, zid] + lost_surface
                            }
                        }

                        # Calculate weighted angles per voxel
                        if (AverageWeightedAngles == 1) {
                            # Calculate angle for the shoot relative to the plane x and y
                            position1 <- c(ShootsBegin$xbegin[j], ShootsBegin$ybegin[j], ShootsBegin$zbegin[j])
                            position2 <- c(ShootsBegin$xend[j], ShootsBegin$yend[j], ShootsBegin$zend[j])
                            V <- position2 - position1
                            alpha <- sum(C*V) / (sqrt(V[1]^2 + V[2]^2 + V[3]^2) * sqrt(C[1]^2 + C[2]^2 + C[3]^2))
                            ShootAngle <- abs(90 - (acos(alpha) / pi*180))  # Angle for shoot

                            # Calculate weighted angle for the voxel
                            tmp1 <- (Mat_surface_per_cell[xid, yid, zid] - ((ShootsBegin$length[j] / (numX*numY*numZ)) * ShootsBegin$diameter[j] * pi/2)) / Mat_surface_per_cell[xid, yid, zid] * Mat_weighted_angle_per_cell[xid, yid, zid]
                            tmp2 <- (((ShootsBegin$length[j] / (numX*numY*numZ)) * ShootsBegin$diameter[j] * pi/2)) / Mat_surface_per_cell[xid, yid, zid] * ShootAngle
                            Mat_weighted_angle_per_cell[xid, yid, zid] <- tmp1 + tmp2
                        }

                    }
                }
            }
        }


        # Loop through all trunks and calculate total surface area, surface loss and weighted angles
        # Get all trees that die during time step
        DeadSegments <- TrunksBegin$treeID[!is.element(TrunksBegin$treeID, TrunksEnd$treeID)]
        DeadSegments1 <- is.element(DeadSegments, TrunksBegin$treeID)
        locDeadSegments <- match(DeadSegments, TrunksBegin$treeID)  # nomatch=0?
        CounterDead <- 1
        TotalDead <- length(locDeadSegments)

        for (j in 1:nrow(TrunksBegin)) {
            X <- ceiling(TrunksBegin$x[j])  # X voxel of tree
            Y <- ceiling(TrunksBegin$y[j])  # Y voxel of tree
            Height <- TrunksBegin$height[j]  # Height of tree
            Diameter <- TrunksBegin$diameter[j]  # Diameter of tree

            SurfaceAreaTotal <- 0

            for (Z in rev(1:ceiling(Height))) {

                # Calculate new total surface area per voxel (for trunks)
                if (TotalSurfaceAreaOpt == 1) {
                    hCone <- Height - Z + 1  # height of cylinder from top to bottom of voxel
                    # rCone <- (hCone / Height) * (Diameter / 2)  # radius of cylinder at bottom of voxel
                    rCone <- Diameter / 2  # radius of cylinder at bottom of voxel

                    # Calculate total surface area in voxel and save it in matrix
                    SurfaceAreaInVoxel <- pi * rCone * sqrt((rCone^2) + (hCone^2)) - SurfaceAreaTotal
                    Mat_surface_per_cell[X, Y, Z] <- Mat_surface_per_cell[X, Y, Z] + SurfaceAreaInVoxel

                    # Update total surface area of cylinder so far (to use in next step)
                    SurfaceAreaTotal <- SurfaceAreaInVoxel + SurfaceAreaTotal
                }

                # If trunk is lost during this time step, add it to lost surface
                if (SurfaceAreaLossOpt == 1) {
                    if (j == locDeadSegments[CounterDead]) {
                        CounterDead <- min(TotalDead, CounterDead + 1)
                        Mat_surfaceloss_per_cell[X, Y, Z] <- Mat_surfaceloss_per_cell[X, Y, Z] + SurfaceAreaInVoxel
                    }
                }

                # Update weighted angle for the voxel
                if (AverageWeightedAngles == 1) {
                    tmp1 <- (Mat_surface_per_cell[X, Y, Z] - SurfaceAreaInVoxel) / Mat_surface_per_cell[X, Y, Z] * Mat_weighted_angle_per_cell[X, Y, Z]
                    tmp2 <- SurfaceAreaInVoxel / Mat_surface_per_cell[X, Y, Z] * 90  # upright 90 degrees angle assumed
                    Mat_weighted_angle_per_cell[X, Y, Z] <- tmp1 + tmp2
                }
            }
        }


        # Calculate light conditions in voxels (relative light conditions)
        if (LightConditionsOpt == 1) {

            # Load file containing information about leaf area per voxel
            voxelsFileName <- paste(voxelFile, i, ".txt", sep="")
            # Voxels <- read.table(file.path(src_dir, voxelsFileName), sep="\t", header=TRUE, skip=1)  # Use this if we fix the trailing tab
            colnames <- as.character(read.table(file.path(src_dir, voxelsFileName), sep="\t", skip=1, nrows=1))
            Voxels <- read.table(file.path(src_dir, voxelsFileName), sep="\t", header=FALSE, skip=2, col.names=append(colnames, "empty_column"))
            Voxels <- Voxels[colnames]

            # Voxel file start with x=y=z=0 => synchronize with matrices used here
            Voxels$x <- Voxels$x + 1
            Voxels$y <- Voxels$y + 1
            Voxels$z <- Voxels$z + 1

            # Store information on leaf area in matrix
            for (j in 1:nrow(Voxels)) {
                Mat_leafArea_per_cell[Voxels$x[j], Voxels$y[j], Voxels$z[j]] <- Voxels$leafarea[j]
            }

            # Calculate single column light conditions based on leaf area distribution
            for (x in 1:dimX) {
                for (y in 1:dimY) {
                    for (z in 1:dimZ) {
                        Mat_light_per_cell[x, y, z] <- exp(-kL * (sum(Mat_leafArea_per_cell[x, y, z:dimZ]) / 10000))
                    }
                }
            }

            # Copy light conditions
            Mat_light_per_cell_Copy <- Mat_light_per_cell

            # Calculate final light conditions by accounting for the light
            # conditions in adjacent voxels
            for (x in corridor:(dimX-corridor)) {
                for (y in corridor:(dimY-corridor)) {
                    for (z in 1:dimZ) {
                        TotalContribution <- 0

                        # loop over ring surrounding the focal voxel
                        for (xx in (x-DistVoxToConsider):(x+DistVoxToConsider)) {
                            for (yy in (y-DistVoxToConsider):(y+DistVoxToConsider)) {
                                Ring <- max(abs(xx - x), abs(yy - y))
                                Contribution <- (1 / (DistVoxToConsider + 1)) * (1 / max(1, (Ring * 8))) * Mat_light_per_cell_Copy[xx, yy, z]
                                TotalContribution <- TotalContribution + Contribution
                            }
                        }

                        Mat_light_per_cell[x, y, z] <- TotalContribution

                    }
                }
            }

        }


        # Store information in Microhabitat matrix and save matrix for this
        # timestep
        # possibly only 5 dimensions to save space
        Microhabitat <- array(rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3] * MatrixDimension), dim=c(dimPlot[1], dimPlot[2], dimPlot[3], MatrixDimension))

        if (TotalSurfaceAreaOpt == 1) {
            Microhabitat[ , , , 1] <- Mat_surface_per_cell[(corridor+1):(dimX-corridor), (corridor+1):(dimY-corridor), 1:dimZ]
        }

        if (SurfaceAreaLossOpt == 1) {
            if (MicrohabitatType == 1) {
                Microhabitat[ , , , 2] <- Mat_surfaceloss_per_cell[(corridor+1):(dimX-corridor), (corridor+1):(dimY-corridor), 1:dimZ] / Mat_surface_per_cell[(corridor+1):(dimX-corridor), (corridor+1):(dimY-corridor), 1:dimZ]
            }

            if (MicrohabitatType == 2) {
                Microhabitat[ , , , 2] <- 0
            }
        }

        if (LightConditionsOpt == 1) {
            Microhabitat[ , , , 3] <- Mat_light_per_cell[(corridor+1):(dimX-corridor), (corridor+1):(dimY-corridor), 1:dimZ]
        }

        if (AverageWeightedAngles == 1) {
            Microhabitat[ , , , 4] <- Mat_weighted_angle_per_cell[(corridor+1):(dimX-corridor), (corridor+1):(dimY-corridor), 1:dimZ]
        }

        MicrohabitatMatSave <- paste("MicrohabitatMatrix", i, ".rds", sep="")
        saveRDS(Microhabitat, file.path(DirectoryMatrices, MicrohabitatMatSave))

        end_time <- Sys.time()
        print(end_time - start_time)

        # flush.console()
        # stop("stop")
    }

    # Save dimensions of plot in seperate file
    saveRDS(dimPlot, file.path(DirectoryMatrices, "dimPlot.rds"))
}
