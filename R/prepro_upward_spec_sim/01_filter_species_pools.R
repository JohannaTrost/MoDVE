library(readr)
library(dplyr)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")
shift_species_q75 <- read_csv(file.path(base_dir, "a5_species_shift_upwards_cc_vs_no_cc.csv"))

base_dir_new <- file.path(base_dir, "a2_1_upward_shifted")
dir.create(base_dir_new, showWarnings = FALSE)

for (pool in 1:10) {
  cat("Processing Species Pool", pool, "\n")
  sp <- read_csv(file.path(base_dir, "a2_1", paste0("SpeciesPool", pool, ".csv")))

  sp_filtered <- shift_species_q75 %>%
    filter(SpeciesPool == pool) %>%
    select(SpeciesID) %>%
    distinct() %>%
    inner_join(., sp, by = "SpeciesID")

  write_csv(sp_filtered, file.path(base_dir_new, paste0("SpeciesPool", pool, ".csv")))

  cat("Filtered species count:", nrow(sp_filtered), "\n")
}

# Copy file to new directory
file.copy(file.path(base_dir, "a2_1", "TraitRanges.csv"),
          file.path(base_dir_new, "TraitRanges.csv"),
          overwrite = TRUE)
