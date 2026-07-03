#' Generate microhabitat matrices from MoF3D output
#'
#' This script integrates forest structure and light conditions from MoF3D
#' into time-dependent microhabitat matrices.
#' It is intended to be executed as a standalone script.
#'
#' @details
#' Usage: Rscript model_pipeline/01_generate_microhabitat.R --config config.toml
#'
#' Example config.toml
#'
#' microhabitatType = 1                            # {1, 2} 1: dynamic forest, 2: static forest (forest at timeStepStart)
#' kL = 0.6                                        # light extinction coefficient
#' DistVoxToConsider = 8                           # no. rings around voxel for light model
#' TotalSurfaceAreaOpt = 1                         # {0, 1} include total surface area
#' SurfaceAreaLossOpt = 1                          # {0, 1} include surface area loss
#' LightNicheOpt = 1                               # {0, 1} include light condidions
#' AverageWeightedAngles = 0                       # {0, 1} include weighted angle averages
#' DirectoryGroIMP = "/path/to/groimp/"            # with MoF3D output "Results/" and "Model/"
#' DirectorySaveMain = "/path/to/output/directory" # e.g. "../modve_data/modve_output/microhabitat"
#' ReplicateForest = 0                             # forest replicate {0, 1, 2}
#' timeStepStart = 1                               # first simulation timestep
#' timeStepEnd = 3                                 # last simulation timestep
#'
#' Output:
#' [microhabitatMatrix<timeStepStart>.rds, ..., microhabitatMatrix<timeStepEnd>.rds]
#' With matrices xDim X yDim X zDim X nVariables,
#' Forest parameter files (Forest_param_global.txt, Forest_param_pass0.txt, dimPlot.rds)
#' are copied to the output directory.
#'
NULL

source("utils.R")

library(data.table)

config <- parse_config()

# ------------------- Parameters ----------------------- #
# Parameters that need to be specified/checked before running this script

# This parameter determines which type of microhatiat matrices are generated:
# 1: real GroIMP forest with dynamics
# 2: static GroIMP forest (only forest at timeStepStart is used)
microhabitatType <- config$microhabitatType

# Parameters of light model
kL <- config$kL  # light extinction coefficient
DistVoxToConsider <- config$DistVoxToConsider  # How many ring around focal voxel to consider in light model (5 voxels in x and y direction)

# Choose the forest parameters that shall be calculated and stored in the microhabitat matrix
# (this list can be extended for possible new applications of the epiphyte model.
# 1: use this variable
# 0: do not use it
TotalSurfaceAreaOpt <- config$TotalSurfaceAreaOpt
SurfaceAreaLossOpt <- config$SurfaceAreaLossOpt
LightNicheOpt <- config$LightNicheOpt
AverageWeightedAngles <- config$AverageWeightedAngles

# Parameters that need to be specified if microhabitat_type=1 or microhabitat_type=2
# Directory of GroIMP files (this directory is stored in the microhabitat folder so that the
# connection to the input GroIMP files is always clear)
DirectoryGroIMP <- config$DirectoryGroIMP
# Directory to save results
DirectorySaveMain <- config$DirectorySaveMain

ReplicateForest <- config$ReplicateForest

# start and end timestep
timeStepStart <- config$timeStepStart
timeStepEnd <- config$timeStepEnd

# Parameters that need to be specified if microhabitat_type=3
# The following parameters are only needed if microhabitat_type=3
# Dimensions of the theoretical forest
# ForestHeight <- 40
# dimXTheoretical <- 50
# dimyTheoretical <- 50
# BAI <- 3  # branch area index for the static, theoretical forest
# LAI <- 6  # leaf are index

# --------------------------------------------------------------------------- #

# Additional parameters
# Names of essential GroIMP files
shootFile <- paste0("shoots_replicate_", ReplicateForest, "_time_step_")
trunkFile <- paste0("trees_replicate_", ReplicateForest, "_time_step_")
voxelFile <- paste0("voxel_replicate_", ReplicateForest, "_time_step_")

C <- c(0, 0, 1)  # Vector orthogonal to plane of X and Y
TotalVoxels <- (DistVoxToConsider * 2 + 1) ^ 2  # Total number of adjacent voxels considered
MatrixDimension <- sum(c(TotalSurfaceAreaOpt, SurfaceAreaLossOpt, LightNicheOpt, AverageWeightedAngles))

# Other parameters created during porting to R
# Those are needed because of GroIMPs output dir structure
model_dir_name <- "Model"
results_dir_name <- "Results"
forest_global_param_name <- "Forest_param_global.txt"
forest_pass_param_name <- paste0("Forest_param_pass", ReplicateForest, ".txt")


# The following parameters are generated automatically
# Load dimensions of forest patch from global forest file
if (microhabitatType == 1 || microhabitatType == 2) {
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
DirectoryMatrices <- file.path(DirectorySaveMain)
dir.create(DirectoryMatrices, recursive=TRUE)

# Copy global and pass forest file to microhabitat folder
if (microhabitatType == 1 || microhabitatType == 2) {
    file.copy(file.path(DirectoryGroIMP, model_dir_name, forest_global_param_name), DirectoryMatrices, overwrite=FALSE)
    file.copy(file.path(DirectoryGroIMP, model_dir_name, forest_pass_param_name), DirectoryMatrices, overwrite=FALSE)
}


# Generation of microhabitat matrix of static or dynamic forest (microhabitat_type=1 or microhabitat_type=2)
# Here, the choosen parameter (total surface, surface loss, light conditions,average angle) are calculated for each voxel in each timestep

if (microhabitatType == 1 || microhabitatType == 2) {

    # In a static forest, only the initial forest at time step timeStepStart is of interest
    if (microhabitatType == 2) {
        timeStepEnd <- timeStepStart
    }

    for (i in int_seq(from=timeStepStart, to=timeStepEnd, by=1)) {

        print(paste("Time step", i))
        start_time <- Sys.time()

        # Load shoot and trunk files of actual and next timestep: Shoots at begin of year
        # and at the end of year/begin of next year
        src_dir <- file.path(DirectoryGroIMP, results_dir_name)

        shootFileNameOld <- paste(shootFile, i, ".txt", sep="")
        shootFileNameNew <- paste(shootFile, i+1, ".txt", sep="")

        trunkFileNameOld <- paste(trunkFile, i, ".txt", sep="")
        trunkFileNameNew <- paste(trunkFile, i+1, ".txt", sep="")

        # In subsequent iterations, we reuse the previous "end data" as the new "begin data"
        # to avoid redundant file loading from disk.
        if (i != timeStepStart) {
            ShootsBegin <- ShootsEnd
            ShootsEnd <- read.table(file.path(src_dir, shootFileNameNew), sep="\t", header=TRUE, skip=1)

            TrunksBegin <- TrunksEnd
            TrunksEnd <- read.table(file.path(src_dir, trunkFileNameNew), sep="\t", header=TRUE, skip=8)
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
        DeadSegments <- ShootsBegin$shootID[!is.element(ShootsBegin$treeID, ShootsEnd$treeID)]
        DeadShootSet <- data.table(shootID=DeadSegments)
        ShootsDT <- data.table(shootID=ShootsBegin$shootID)
        ShootsDT <- ShootsDT[, is_dead := shootID %in% DeadShootSet$shootID] # faster than match
        locDeadSegments <- which(ShootsDT$is_dead)
        CounterDead <- 1
        TotalDead <- length(locDeadSegments)

        # Loop through all shoots and calculate total surface area, surface loss and weighted angles
        for (j in seq_len(nrow(ShootsBegin))) {

            # Get voxel the shoot is intersecting with
            # Here, we are using a simple approach because the shoots are
            # usually small only intersect with maximum two voxels
            UniqueX <- c(ceiling(ShootsBegin$xbegin[j]), ceiling(ShootsBegin$xend[j]))
            UniqueY <- c(ceiling(ShootsBegin$ybegin[j]), ceiling(ShootsBegin$yend[j]))
            UniqueZ <- c(ceiling(ShootsBegin$zbegin[j]), ceiling(ShootsBegin$zend[j]))

            numX <- sum((UniqueX[2] - UniqueX[1]) > 0, na.rm=TRUE) + 1
            numY <- sum((UniqueY[2] - UniqueY[1]) > 0, na.rm=TRUE) + 1
            numZ <- sum((UniqueZ[2] - UniqueZ[1]) > 0, na.rm=TRUE) + 1

            for (x in seq_len(numX)) {
                xid <- UniqueX[x]

                for (y in seq_len(numY)) {
                    yid <- UniqueY[y]

                    for (z in seq_len(numZ)) {
                        zid <- UniqueZ[z]

                        # Calculate new total surface area per voxel
                        if (TotalSurfaceAreaOpt == 1) {
                            # Surface of single branch within voxel
                            voxel_branch_surface <- (ShootsBegin$length[j] / (numX*numY*numZ)) * ShootsBegin$diameter[j] * pi/2

                            Mat_surface_per_cell[xid, yid, zid] <- Mat_surface_per_cell[xid, yid, zid] + voxel_branch_surface
                        }

                        # If branch is lost during this time step, add it to lost surface
                        if (SurfaceAreaLossOpt == 1) {
                            if (j == locDeadSegments[CounterDead] & TotalDead > 0) {
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
        DeadTrunkSet <- data.table(treeID=DeadSegments)
        TrunksDT <- data.table(treeID=TrunksBegin$treeID)
        TrunksDT <- TrunksDT[, is_dead := treeID %in% DeadTrunkSet$treeID] # faster than match
        locDeadSegments <- which(TrunksDT$is_dead)
        CounterDead <- 1
        TotalDead <- length(locDeadSegments)

        for (j in seq_len(nrow(TrunksBegin))) {
            X <- ceiling(TrunksBegin$x[j])  # X voxel of tree
            Y <- ceiling(TrunksBegin$y[j])  # Y voxel of tree
            Height <- TrunksBegin$height[j]  # Height of tree
            Diameter <- TrunksBegin$diameter[j]  # Diameter of tree

            SurfaceAreaTotal <- 0

            for (Z in rev(seq_len(ceiling(Height)))) {

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
                print(paste0("Total dead: ", TotalDead))
                if (SurfaceAreaLossOpt == 1 & TotalDead > 0) {
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
        if (LightNicheOpt == 1) {

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
            Mat_leafArea_per_cell[cbind(Voxels$x, Voxels$y, Voxels$z)] <- Voxels$leafarea

            # Calculate single column light conditions based on leaf area distribution
            for (x in seq_len(dimX)) {
                for (y in seq_len(dimY)) {
                    for (z in seq_len(dimZ)) {
                        Mat_light_per_cell[x, y, z] <- exp(-kL * (sum(Mat_leafArea_per_cell[x, y, z:dimZ]) / 10000))
                    }
                }
            }

            # Copy light conditions
            Mat_light_per_cell_Copy <- Mat_light_per_cell

            # Calculate final light conditions by accounting for the light conditions in adjacent voxels
            for (x in int_seq(from=corridor, to=dimX-corridor, by=1)) {
                for (y in int_seq(from=corridor, to=dimY-corridor, by=1)) {
                    for (z in seq_len(dimZ)) {
                        TotalContribution <- 0

                        # loop over ring surrounding the focal voxel
                        for (xx in int_seq(from=x-DistVoxToConsider, to=x+DistVoxToConsider, by=1)) {
                            for (yy in int_seq(from=y-DistVoxToConsider, to=y+DistVoxToConsider, by=1)) {
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

        # Store information in microhabitat matrix and save matrix for this timestep
        microhabitat <- array(rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3] * (MatrixDimension+1)),
                              dim=c(dimPlot[1], dimPlot[2], dimPlot[3], (MatrixDimension+1)))

        idx1 <- int_seq(from=corridor+1, to=dimX-corridor, by=1)
        idx2 <- int_seq(from=corridor+1, to=dimY-corridor, by=1)
        idx3 <- seq_len(dimZ)

        if (TotalSurfaceAreaOpt == 1) {
            # Compute PAI
            Mat_plantArea_per_cell <- Mat_surface_per_cell + (Mat_leafArea_per_cell / 10000)
            microhabitat[ , , , MatrixDimension+1] <- Mat_plantArea_per_cell[idx1, idx2, idx3]
            microhabitat[ , , , 1] <- Mat_surface_per_cell[idx1, idx2, idx3]
        }

        if (SurfaceAreaLossOpt == 1) {
            if (microhabitatType == 1) {
                microhabitat[ , , , 2] <- Mat_surfaceloss_per_cell[idx1, idx2, idx3] / Mat_surface_per_cell[idx1, idx2, idx3]
            }

            if (microhabitatType == 2) {
                microhabitat[ , , , 2] <- 0
            }
        }

        if (LightNicheOpt == 1) {
            microhabitat[ , , , 3] <- Mat_light_per_cell[idx1, idx2, idx3]
        }

        if (AverageWeightedAngles == 1) {
            microhabitat[ , , , 4] <- Mat_weighted_angle_per_cell[idx1, idx2, idx3]
        }

        microhabitatMatSave <- paste("microhabitatMatrix", i, ".rds", sep="")
        saveRDS(microhabitat, file.path(DirectoryMatrices, microhabitatMatSave))

        end_time <- Sys.time()
        print(end_time - start_time)
    }

    # Save dimensions of plot in seperate file
    saveRDS(dimPlot, file.path(DirectoryMatrices, "dimPlot.rds"))
}
