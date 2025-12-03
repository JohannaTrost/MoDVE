config <- parse_config("tests/config_a1.toml")

# only for MicrohabitatType 1 or 2

read_vox_dt <- function(voxelsFileName) {

  colnames <- voxelsFileName |>
    utils::read.table(sep = "\t", skip = 1, nrows = 1) |>
    as.character()

  Voxels <- voxelsFileName |>
    read.table(sep = "\t", header = FALSE, skip = 2,
               col.names = append(colnames, "empty_column"))
  Voxels <- Voxels[colnames]

  # Voxel file start with x=y=z=0 => synchronize with matrices used here
  # TODO: this is quite bad, MoF3D outputs need to be made consistent
  Voxels$x <- Voxels$x + 1
  Voxels$y <- Voxels$y + 1
  Voxels$z <- Voxels$z + 1

  return(Voxels)
}

{
  # ------------------- Parameters ----------------------- #
  # Parameters that need to be specified/checked before running this script

  # This parameter determines which type of microhatiat matrices are generated:
  # 1: real GroIMP forest with dynamics
  # 2: static GroIMP forest (only forest at timeStepStart is used)
  MicrohabitatType <- config$MicrohabitatType

  # Parameters that need to be specified if MicrohabitatType=1 or MicrohabitatType=2
  # Directory of GroIMP files (this directory is stored in the Microhabitat folder so that the
  # connection to the input GroIMP files is always clear)
  mof3d_dir <- config$DirectoryGroIMP
  # Directory to save results
  DirectorySaveMain <- config$DirectorySaveMain

  rep_forest <- config$ReplicateForest

  # start and end timestep
  timeStepStart <- config$timeStepStart
  timeStepEnd <- config$timeStepEnd
}

# Other parameters created during porting to R
# Those are needed because of GroIMPs output dir structure
mof3d_res_dir <- file.path(mof3d_dir, "Results")
mof3d_model_dir <- file.path(mof3d_dir, "Model")
forest_global_file <- "Forest_param_global.txt"
forest_pass_file <- paste("Forest_param_pass", rep_forest, ".txt", sep = "")

# Create folder to save the microhabitat matrices
output_dir <- file.path(DirectorySaveMain)
dir.create(output_dir, recursive = TRUE)
# Copy Mof3D forest data to result directory
file.path(mof3d_model_dir, forest_global_file) |>
  file.copy(to = output_dir, overwrite = FALSE)
file.path(mof3d_model_dir, forest_pass_file) |>
  file.copy(to = output_dir, overwrite = FALSE)

# Load dimensions of forest patch from global forest file
path_to_forest_global <- file.path(mof3d_model_dir, forest_global_file)
GlobalForest <- read.table(path_to_forest_global, sep = "\t", row.names = 1)
config$MaxX <- GlobalForest["MaxX", 1]
config$MaxY <- GlobalForest["MaxY", 1]
config$MaxZ <- GlobalForest["MaxZ", 1]
config$corridor <- GlobalForest["WidthCorridor", 1]

shootFile <- paste("shoots_replicate_", rep_forest, "_time_step_", sep = "")
trunkFile <- paste("trees_replicate_", rep_forest, "_time_step_", sep = "")
voxFile <- paste("voxel_replicate_", rep_forest, "_time_step_", sep = "")

# Load data for first time step
ShootsBegin <- file.path(mof3d_res_dir, paste(shootFile, timeStepStart, ".txt", sep = "")) |>
  read.table(sep = "\t", header = TRUE, skip = 1)
TrunksBegin <-   file.path(mof3d_res_dir, paste(trunkFile, timeStepStart, ".txt", sep = "")) |>
  read.table(sep = "\t", header = TRUE, skip = 8)

if (MicrohabitatType == 2) timeStepEnd <- timeStepStart

for (i in timeStepStart:timeStepEnd) {

  print(paste("Time step", i))

  # Load shoot and trunk files of actual and next timestep: Shoots at begin of year
  # and at the end of year/begin of next year
  # Previous time step
  if (i != timeStepStart) {
    # Last step end becomes beginning of this step
    ShootsBegin <- ShootsEnd
    TrunksBegin <- TrunksEnd
  }
  ShootsEnd <- paste(shootFile, i + 1, ".txt", sep = "") |>
    read.table(sep = "\t", header = TRUE, skip = 1)
  TrunksEnd <- paste(trunkFile, i + 1, ".txt", sep = "") |>
    read.table(sep = "\t", header = TRUE, skip = 8)

  # Find which branch segments and trunks die this time step
  dead_branches_id <- find_dead_segments(ShootsBegin, ShootsEnd)
  dead_trees_id <- find_dead_trees(TrunksBegin, TrunksEnd)

  # Load file containing information about leaf area per voxel
  if (config$LightConditionsOpt) {
    vox_dt <- file.path(mof3d_res_dir, paste(voxFile, i, ".txt", sep = "")) |>
      read_vox_dt()
  }

  path_to_output <- file.path(
    output_dir, paste("MicrohabitatMatrix", i, ".rds", sep = "")
  )

  create_microhabitat_mat(
    config = config,
    shoot_dt = ShootsBegin,
    trunk_dt = TrunksBegin,
    vox_dt = ifelse(config$LightConditionsOpt, vox_dt, NULL),
    path_to_output = path_to_output,
    dead_branches_id = dead_branches_id,
    dead_trees_id = dead_trees_id
    )

}
