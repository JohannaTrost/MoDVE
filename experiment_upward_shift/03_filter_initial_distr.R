# -----
# Filter initial spatial distribution keeping only species that shift upward under CC

library(readr)
library(dplyr)

base_dir <- file.path("../modve_data/modve_output/regua")
zenodo_base_dir <- file.path("../modve_data_zenodo/modve_output/regua")
shift_species_q75 <- read_csv(file.path(base_dir, "species_shift_upwards_cc_vs_no_cc.csv"))

base_dir_new <- file.path(base_dir, "climdata_era5_cmip6_1981-2100_ssp245_no_cc", "upward_shifted_distribution")
dir.create(base_dir_new, showWarnings = FALSE)

for (forest in 0:2) {
  dir.create(file.path(base_dir_new, paste0("forest", forest)), showWarnings = FALSE)

  for (pool in 1:10) {
    path <- file.path(zenodo_base_dir, "climdata_era5_cmip6_1981-2100_ssp245_no_cc", "distribution",
                      paste0("forest", forest), paste0("ID_SpeciesP_", pool, "_Rep_1.csv"))
    initialDistr <- read_csv(path, show_col_types = FALSE)

    sp_filtered <- shift_species_q75 %>%
        filter(SpeciesPool == pool) %>%
        select(SpeciesID) %>%
        distinct() %>%
        inner_join(., initialDistr, by = "SpeciesID") %>%
        arrange(IndividualID) %>%
        mutate(IndividualID = 1:n())

    out_path <- file.path(base_dir_new, paste0("forest", forest),
                          paste0("ID_SpeciesP_", pool, "_Rep_1.csv"))
    write_csv(sp_filtered, out_path)
  }
}



