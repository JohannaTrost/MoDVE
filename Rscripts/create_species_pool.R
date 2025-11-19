
species_params <- parse_config("tests/config_a2.toml")
check_species_params(species_params)

# Prepare output table
SaveDirectory <- file.path(species_params$MainOutputDirectory)
dir.create(SaveDirectory, recursive=TRUE)
trait_names <- species_trait_names()
sp_traits_df <- data.frame(matrix(ncol = length(trait_names), nrow = 0))
colnames(sp_traits_df) <- trait_names
# TODO: careful, column types are all logical instead of numeric

for (pool_nb in 1:species_params$numSpeciesPools) {
  for (sp in 1:species_params$NumberOfSpecies) {

    species_traits <- draw_species_traits(species_params)
    species_traits$SpeciesID <- sp
    sp_traits_df[sp, ] <- as.data.frame(species_traits)
  }

  # Save trait dataframe
  SpeciesPoolFileName <- paste("SpeciesPool", pool_nb, ".csv", sep = "")
  write.csv(sp_traits_df, file.path(SaveDirectory, SpeciesPoolFileName),
            row.names = FALSE)
}
