# -----
# Species pools with upward shifting species under CC

library(readr)
library(dplyr)

base_dir <- file.path("../modve_data/modve_output/regua")
zenodo_base_dir <- file.path("../modve_data_zenodo/modve_output/regua")
shift_species_q75 <- read_csv(file.path(base_dir, "species_shift_upwards_cc_vs_no_cc.csv"))

base_dir_new <- file.path(base_dir, "upward_shifted_species_pool")
dir.create(base_dir_new, showWarnings = FALSE)

for (pool in 1:10) {
  cat("Processing Species Pool", pool, "\n")
  sp <- read_csv(file.path(zenodo_base_dir, "species_pools", paste0("SpeciesPool", pool, ".csv")))

  sp_filtered <- shift_species_q75 %>%
    filter(SpeciesPool == pool) %>%
    select(SpeciesID) %>%
    distinct() %>%
    inner_join(., sp, by = "SpeciesID")

  write_csv(sp_filtered, file.path(base_dir_new, paste0("SpeciesPool", pool, ".csv")))

  cat("Filtered species count:", nrow(sp_filtered), "\n")
}

# Copy file to new directory TODO update path
file.copy(file.path(zenodo_base_dir, "species_pools", "TraitRanges.csv"),
          file.path(base_dir_new, "TraitRanges.csv"),
          overwrite = TRUE)
