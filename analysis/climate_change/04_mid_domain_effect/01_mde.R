# -----
# Run mid-domain effect simulation
# 1. Extract simulated species ranges
# 2. Randomly place ranges along the vertical profile 1000 times
# 3. Get species vertical species richness distribution

library(readr)
library(dplyr)
library(purrr)

# Set year for MDE simulation
year <- 2050 # Do for different years

# --- Mid domain randomization to compare to actual richness

# Load ranges
base_dir <- "../modve_data/modve_output/regua"
spec_distr <- read_csv(file.path(base_dir, "species_distribution_cc_vs_no_cc.csv"))

# --- Step 1. Extract realized vertical ranges and max height in forest ---

spec_distr_stats <- spec_distr %>%
  filter(Scenario == "CC", Year == year) %>%
  group_by(SpeciesPool, SpeciesID, ForestID) %>%
  summarize(Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop")

# Compute max height per SpeciesPool × ForestID
max_height_df <- spec_distr %>%
  filter(Scenario == "CC") %>%
  group_by(SpeciesPool, ForestID) %>%
  summarize(MaxHeight = max(Height, na.rm = TRUE), .groups = "drop")

# Add it to spec_distr_stats
spec_distr_stats <- spec_distr_stats %>%
  left_join(max_height_df, by = c("SpeciesPool", "ForestID"))

# --- Step 2. Function to run one MDR ---
run_mdr <- function(speciesPool, forest) {
  # get scalar max height for this pool × forest
  max_h_scalar <- max_height_df %>%
    filter(SpeciesPool == speciesPool, ForestID == forest) %>%
    pull(MaxHeight)

  bins <- 0:ceiling(max_h_scalar)

  # Filter to the species in this pool × forest
  vr <- spec_distr_stats %>%
    filter(SpeciesPool == speciesPool, ForestID == forest) %>%
    mutate(
      MaxHeightSpec = MaxHeight - Range
    ) %>%
    # vectorized safe random placement: uniform * MaxHeightSpec
    mutate(
      min_h = runif(n()) * MaxHeightSpec,
      max_h = min_h + Range
    )

  map_dfr(seq_along(bins[-length(bins)]), function(i) {
    bin_min <- bins[i]
    bin_max <- bins[i+1]

    n_species <- vr %>%
      filter(min_h < bin_max, max_h > bin_min) %>%
      distinct(SpeciesID) %>%            # count distinct species
      nrow()

    tibble(bin_mid = (bin_min + bin_max)/2,
           Richness = n_species)
  })
}

# --- Step 3. Run many randomizations (Monte Carlo) ---

set.seed(42)
n_iter <- 1000

all_results <- map_dfr(1:10, function(sp) {
  map_dfr(0:2, function(f) {
    map_dfr(1:n_iter, ~ run_mdr(sp, f) %>%
              mutate(iter = .x)) %>%
      mutate(ForestID = f)
  }) %>%
    mutate(SpeciesPool = sp)
})

# --- Step 4. Aggregate across randomizations ---

expected_curve <- all_results %>%
  group_by(bin_mid) %>%
  summarise(
    mean_richness = mean(Richness, na.rm = TRUE),
    sd_richness   = sd(Richness, na.rm = TRUE),
    n_obs         = n(),
    sem_richness  = sd_richness / sqrt(n_obs),
    .groups = "drop"
  )
write_csv(expected_curve,
          file.path(base_dir, "climdata_era5_cmip6_1981-2100_ssp245", paste0("mde_richness_", year, ".csv")))