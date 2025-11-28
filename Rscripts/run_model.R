options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

# Parse input configuration file
config <- parse_config("tests/config_a4.toml")

# Output directory
dir_output <- config$DirectoryModelResults
dir.create(dir_output, recursive = TRUE)

# Detect the number of CPU cores and register the parallel backend
# Note: detectCores() will detect the total number of cores on a HPC node,
# so use the $SLURM_NTASKS env var if defined.
nTasks <- Sys.getenv("SLURM_NTASKS")
if (nTasks != "") {
  numCores <- strtoi(nTasks)
} else {
  numCores <- parallel::detectCores()
}
doParallel::registerDoParallel(numCores)

# RNG seed
random_state_file <- config$RandomState
if (!is.null(random_state_file) && file.exists(random_state_file)) {
  writeLines("Loading previous random number generator state from ", random_state_file)
  restore_rng(random_state_file)
} else {
  # If we provide an integer, use it to set the seed otherwise it will
  # be NULL and therefore a random seed will be created.
  seed <- config$seed
  set.seed(seed, kind = "Mersenne-Twister")
}
# Save current random state.
writeLines("Storing the random number generator state...")
save_rng(file.path(dir_output, "random_state_seed.RData"))

# Main loop for the community model for each species pool and for each replicate
# Each loop employs a different random number generator (RNG) stream, resulting in distinct,
# statistically independent random sequences. These sequences are reproducible across multiple
# runs, provided the master seed and other input parameters remain unchanged. Importantly, the
# number of cores does not influence the sequences, only the order in which the loops are
# processed. Consequently, parallel and serial execution will yield identical results.
# Internally, the foreach package employs the L'Ecuyer-CMRG RNG algorithm for reliable random
# number generation, ensuring reproducible results even in parallel computing environments.
# Choose species pools to use and number of replicates per species pool
pairs <- expand.grid(config$numSpeciesPools, config$replicatePerSpeciesPool)
colnames(pairs) <- c("numPool", "r")
output <- foreach::foreach(pair_idx = seq_len(nrow(pairs))) %dorng% {

  numPool <- pairs$numPool[pair_idx]
  r <- pairs$r[pair_idx]

  # Load initial distribution of individuals
  init_filename <- paste("ID_SpeciesP_", numPool, "_Rep_", r, ".csv", sep = "")
  path_to_init_file <- file.path(config$DirectoryModelMain, init_filename)

  # Load species pool
  sp_filename <- paste("SpeciesPool", numPool, ".csv", sep = "")
  SpeciesPoolFileName <- file.path(config$DirectorySpeciesPools, sp_filename)

  # Load microhabitat (static forest)
  Microhabitat <- readRDS(file.path(config$DirectoryMicrohabitat, "MicrohabitatMatrix1.rds"))

  # Create Save-Directory for each each replicate/initialDistribution
  dir_output_this_run <- file.path(dir_output, paste("ID_SpeciesP_", numPool, "_Rep_", r, sep=""))
  dir.create(dir_output_this_run, recursive = TRUE)
  path_to_ind_output <- file.path(dir_output_this_run, "IndividualMatrixTimeStep.csv")
  path_to_sp_output <- file.path(dir_output_this_run, "SpeciesSummary.csv")
  path_to_comm_output <- file.path(dir_output_this_run, "CommunitySummary.csv")

  # Run the IBM
  run_modve_sim(
    config,
    SpeciesPoolFilename,
    Microhabitat,
    path_to_init_file,
    path_to_ind_output,
    path_to_sp_output,
    path_to_comm_output
  )
  return(NULL)
}
