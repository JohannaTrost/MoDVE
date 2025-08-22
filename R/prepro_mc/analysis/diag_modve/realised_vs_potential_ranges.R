library(tidyverse)

# Paths
model_path <- file.path("/Users/johanna/Uni/masterarbeit/data/a5_output/v8_real_niches_original_model_light_hum_temp/ID_SpeciesP_1_Rep_1")
res_file <- "IndividualMatrixTimeStep199.csv"

out_dir <- file.path("../../../data/a4_output/v8_real_niches_original_model_light_hum_temp")
suit_file <- "ScaledSuitability_1_TimeStep199.h5"

 # --- Realised ranges

res <- read_csv(file.path(model_path, res_file))
realised_ranges <- res %>%
  filter(Status == 1) %>%  # Filter for individuals that survived
  mutate(height = Z - 0.5) %>%
  select(SpeciesID, height) %>%
  group_by(SpeciesID) %>%
  summarise(realised_height_range = max(height) - min(height), .groups = "drop")

max_height <- max(res$Z -0.5)

# --- Potential ranges
suitability <- rhdf5::h5read(file.path(out_dir, suit_file), "ScaledSuitabilityScores")

potential_richness <- apply(suitability > quantile(suitability[suitability > 0], 0.1),
                            c(3, 4), sum) > 0
potential_richness_bounded <- potential_richness[1:ceiling(max_height),]

# suppose rows correspond to heights 1:80
heights <- 1:nrow(potential_richness_bounded) - 0.5

# apply over columns
height_ranges <- apply(potential_richness_bounded, 2, function(sp) {
  if (all(sp == 0)) {
    return(NA)  # species never occurs
  } else {
    present_heights <- heights[sp == 1]
    return(max(present_heights) - min(present_heights))
  }
})

potential_ranges_bounded <- data.frame(
  SpeciesID = 1:ncol(potential_richness_bounded),
  potential_height_range = height_ranges
)

# Combine potential and realised ranges

height_ranges <- realised_ranges %>%
    left_join(potential_ranges_bounded, by = "SpeciesID") %>%
    mutate(
        range_filling = realised_height_range / potential_height_range
    ) %>%
    # Sort by range filling
    arrange(desc(range_filling))

# Plot the distribution of all potential ranges
pdf("../../../figs/a4_potential_ranges_bounded.pdf")
ggplot(height_ranges, aes(x = potential_height_range)) +
    geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
    labs(
        x = "Potential Height Range (m)",
        y = "Frequency"
    ) +
    theme_minimal()
dev.off()

# --- MDE potential niches

heights <- 1:nrow(potential_richness) - 0.5

# apply over columns
height_ranges <- apply(potential_richness, 2, function(sp) {
  if (all(sp == 0)) {
    return(NA)  # species never occurs
  } else {
    present_heights <- heights[sp == 1]
    return(max(present_heights) - min(present_heights))
  }
})

potential_ranges <- data.frame(
  SpeciesID = 1:ncol(potential_richness),
  potential_height_range = height_ranges
) %>%
  filter(SpeciesID %in% realised_ranges$SpeciesID)

# --- Compute MDE

max_height <- max(heights)
bins <- 0:ceiling(max_height)

# --- Step 2. Function to run one MDR ---
run_mdr <- function() {
  vr <- potential_ranges %>%
    mutate(max_height = max_height - potential_height_range,
           min_h = runif(n(), min = 0, max = max_height),
           max_h = min_h + potential_height_range)

  map_dfr(seq_along(bins[-length(bins)]), function(i) {
    bin_min <- bins[i]
    bin_max <- bins[i+1]

    n_species <- vr %>%
      filter(min_h < bin_max, max_h > bin_min) %>%
      nrow()

    tibble(bin_mid = (bin_min + bin_max)/2,
           Richness = n_species)
  })
}

# --- Run many randomizations (Monte Carlo) ---
set.seed(42)
n_iter <- 1000

null_results <- map_dfr(1:n_iter, ~run_mdr() %>% mutate(iter = .x))

# --- Step 4. Aggregate across randomizations ---
expected_curve <- null_results %>%
  group_by(bin_mid) %>%
  summarise(
    mean_richness = mean(Richness),
    sd_richness   = sd(Richness),
    .groups = "drop"
  ) %>%
  rename(height = bin_mid)

# Plot the expected curve

speciesRichnessPlot <- ggplot(expected_curve, aes(x = mean_richness, y = height)) +
  geom_path(aes(x = mean_richness, y = height, color = "Avg. MDE Richness"), size = 1, linetype = "dashed") +
  geom_ribbon(aes(y = height, xmin = mean_richness - sd_richness, xmax = mean_richness + sd_richness, fill = "MDE +-SD"), alpha = 0.2) +
  scale_color_manual(name = "Richness Type", values = c("Actual Richness" = "darkblue", "Avg. MDE Richness" = "red")) +
  scale_fill_manual(name = "MDE +-SD", values = c("MDE +-SD" = "lightgrey")) +
  labs(
    x = "Species Richness",
    y = "Height (m)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),       # tick label size
    axis.title = element_text(size = 14)       # axis titles larger
  )

filename <- paste0(DirectoryPlots, "Replicate_1_SpeciesPool_1_MDnull_PotentialVerticalSpeciesRichness_", ts, ".pdf")

pdf(filename)
print(speciesRichnessPlot)
dev.off()