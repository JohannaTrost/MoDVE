
config <- parse_config("tests/config_a3.toml")
set.seed(config$seed, kind="Mersenne-Twister")

# Prepare folder paths
DirectoryModelMain <- config$DirectoryModelMain
DirectoryMicrohabitatMain <- config$DirectoryMicrohabitat
DirectorySpeciesPoolsMain <- config$DirectorySpeciesPools

dir.create(DirectoryModelMain, recursive=TRUE)

# Load microhabitat matrix
microhabitat_filename <- paste("MicrohabitatMatrix", TimeStep, ".rds", sep="")
FileInitalMatrix <- file.path(DirectoryMicrohabitatMain, microhabitat_filename)
Microhabitat <- readRDS(FileInitalMatrix)
# Transform light values from relative to absolute
Microhabitat[, , , 3] <- Microhabitat[, , , 3] * config$Imax

for (numPool in config$numSpeciesPools[1]:config$numSpeciesPools[2]) {

  # Load the species pool
  species_filename <- paste("SpeciesPool", numPool, ".csv", sep = "")
  Input_file <- file.path(DirectorySpeciesPoolsMain, species_filename)
  SpeciesPool <- read.csv(Input_file)

  for (numReplicates in seq_len(config$replicatePerSpeciesPool)) {

    path_to_output <- file.path(DirectoryModelMain, paste(
      "ID_SpeciesP_", numPool, "_Rep_", numReplicates, ".csv", sep = ""
    ))

    # Generate the distribution and individual traits
    draw_initial_individuals(
      config, SpeciesPool, Microhabitat, path_to_output,
      largest_inds_first = config$SingleSpeciesModel,
      largest_voxels_first = config$MethodVoxel
    )
  }
}
