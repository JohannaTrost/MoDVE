library(tidyverse)

# Paths
model_path <- file.path("/Users/johanna/Uni/masterarbeit/data/a5_output/v8_real_niches_original_model_light_hum_temp/ID_SpeciesP_1_Rep_1")
res_file <- "IndividualMatrixTimeStep199.csv"

out_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/a4_output/v8_real_niches_original_model_light_hum_temp")
suit_file <- "ScaledSuitability_1_TimeStep199.h5"

mh_path <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_input/v2_microhabitat_lowland_zach_mh_2024_mc")
mh_file <- "microhabitatMatrix199.rds"

 # --- Realised ranges

res <- read_csv(file.path(model_path, res_file))
realised_ranges <- res %>%
  filter(Status == 1) %>%  # Filter for individuals that survived
  mutate(height = Z - 0.5) %>%
  select(SpeciesID, height) %>%
  group_by(SpeciesID) %>%
  summarise(realised_height_range = max(height) - min(height), .groups = "drop")

max_veg_height <- max(res$Z -0.5)

# --- Potential ranges
suitability <- rhdf5::h5read(file.path(out_dir, suit_file), "ScaledSuitabilityScores")
max_height <- dim(suitability)[3]

potential_richness <- apply(suitability > quantile(suitability[suitability > 0], 0.1),
                            c(3, 4), sum) > 0
full_potential_richness <- potential_richness
potential_richness <- potential_richness[1:ceiling(max_veg_height),]

# suppose rows correspond to heights 1:80
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
)

# Combine potential and realised ranges

height_ranges <- realised_ranges %>%
    left_join(potential_ranges, by = "SpeciesID") %>%
    mutate(
        range_filling = realised_height_range / potential_height_range
    ) %>%
    # Sort by range filling
    arrange(desc(range_filling))


# Reshape to long format for plotting
height_ranges_long <- height_ranges |>
  pivot_longer(cols = c(realised_height_range, potential_height_range),
               names_to = "type", values_to = "height_range")

# Plot the distribution of all potential ranges
pdf("../../figs/a5_plots_test/v8_real_niches_original_model_light_hum_temp/a4_range_distribution.pdf")
ggplot(height_ranges_long, aes(x = height_range, fill = type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Realised and Potential Height Ranges",
       x = "Height Range", y = "Density") +
  theme_minimal()
dev.off()

# --- 1. Compute Bounded MDEs

max_height <- max(heights)
bins <- 0:ceiling(max_height)

# --- Step 1.1 Randomization for Potential ranges ---
run_mdr_potential <- function() {
  vr <- height_ranges %>%
    mutate(
      max_shift = max_height - potential_height_range,
      min_h = runif(n(), min = 0, max = max_shift),
      max_h = min_h + potential_height_range
    )

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

# --- Step 1.2 Randomization for Realised ranges ---
run_mdr_realised <- function() {
  vr <- height_ranges %>%
    mutate(
      max_shift = max_height - realised_height_range,
      min_h = runif(n(), min = 0, max = max_shift),
      max_h = min_h + realised_height_range
    )

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

# --- Step 1.3 Run Monte Carlo for both ---
set.seed(42)
n_iter <- 1000

null_potential <- map_dfr(1:n_iter, ~run_mdr_potential() %>% mutate(iter = .x))
null_realised  <- map_dfr(1:n_iter, ~run_mdr_realised()  %>% mutate(iter = .x))

# --- Step 4. Aggregate into one expected_curve ---
expected_curve <- null_potential %>%
  group_by(bin_mid) %>%
  summarise(
    mean_potential_richness = mean(Richness),
    sd_potential_richness   = sd(Richness),
    .groups = "drop"
  ) %>%
  left_join(
    null_realised %>%
      group_by(bin_mid) %>%
      summarise(
        mean_realised_richness = mean(Richness),
        sd_realised_richness   = sd(Richness),
        .groups = "drop"
      ),
    by = "bin_mid"
  ) %>%
  rename(height = bin_mid)

# --- 1.4 Plot the expected curve ---

speciesRichnessPlot <- ggplot(expected_curve, aes(y = height)) +
  # Potential MDE curve + ribbon
  geom_path(aes(x = mean_potential_richness, color = "Potential MDE"),
            size = 1, linetype = "dashed") +
  geom_ribbon(aes(xmin = mean_potential_richness - sd_potential_richness,
                  xmax = mean_potential_richness + sd_potential_richness,
                  fill = "Potential ±SD"), alpha = 0.2) +

  # Realised MDE curve + ribbon
  geom_path(aes(x = mean_realised_richness, color = "Realised MDE"),
            size = 1, linetype = "dashed") +
  geom_ribbon(aes(xmin = mean_realised_richness - sd_realised_richness,
                  xmax = mean_realised_richness + sd_realised_richness,
                  fill = "Realised ±SD"), alpha = 0.2) +

  scale_color_manual(
    name = "MDE Curve",
    values = c("Potential MDE" = "red", "Realised MDE" = "blue")
  ) +
  scale_fill_manual(
    name = "MDE ±SD",
    values = c("Potential ±SD" = "pink", "Realised ±SD" = "lightblue")
  ) +
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

DirectoryPlots <- file.path("../../figs/a5_plots_test/v8_real_niches_original_model_light_hum_temp")
filename <- file.path(DirectoryPlots,
                      "Replicate_1_SpeciesPool_1_MDnull_Potential_vs_Realised_199.pdf")

pdf(filename)
print(speciesRichnessPlot)
dev.off()


# --- 2. Unbounded MDEs

# --- Step 2.1 Compute full potential ranges ---
heights <- 1:nrow(full_potential_richness) - 0.5

# apply over columns
full_height_ranges <- apply(full_potential_richness, 2, function(sp) {
  if (all(sp == 0)) {
    return(NA)  # species never occurs
  } else {
    present_heights <- heights[sp == 1]
    return(max(present_heights) - min(present_heights))
  }
})

full_potential_ranges <- data.frame(
  SpeciesID = 1:ncol(full_potential_richness),
  potential_height_range = full_height_ranges
)

# Combine potential and realised ranges
full_height_ranges <- realised_ranges %>%
    left_join(full_potential_ranges, by = "SpeciesID") %>%
    mutate(
        range_filling = realised_height_range / potential_height_range
    ) %>%
    # Sort by range filling
    arrange(desc(range_filling))


# --- Step 2.2 Randomization for Potential ranges ---
max_height <- max(heights)
bins <- 0:ceiling(max_height)
run_unbound_mdr_potential <- function() {
  vr <- full_height_ranges %>%
    mutate(
      midpoint = runif(n(), min = 0, max = max_veg_height),
      min_h = midpoint - potential_height_range / 2,
      max_h = min_h + potential_height_range
    )

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

# --- Step 2.3 Randomization for Realised ranges ---
run_mdr_realised <- function() {
  vr <- full_height_ranges %>%
    mutate(
      midpoint = runif(n(), min = 0, max = max_veg_height),
      min_h = midpoint - realised_height_range / 2,
      max_h = min_h + realised_height_range
    )

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

# --- Step 1.3 Run Monte Carlo for both ---
set.seed(42)
n_iter <- 1000

null_potential_unbound <- map_dfr(1:n_iter, ~run_mdr_potential() %>% mutate(iter = .x))
null_realised_unbound  <- map_dfr(1:n_iter, ~run_mdr_realised()  %>% mutate(iter = .x))

# --- Step 4. Aggregate into one expected_curve ---
expected_curve_unbound <- null_potential_unbound %>%
  group_by(bin_mid) %>%
  summarise(
    mean_potential_richness = mean(Richness),
    sd_potential_richness   = sd(Richness),
    .groups = "drop"
  ) %>%
  left_join(
    null_realised_unbound %>%
      group_by(bin_mid) %>%
      summarise(
        mean_realised_richness = mean(Richness),
        sd_realised_richness   = sd(Richness),
        .groups = "drop"
      ),
    by = "bin_mid"
  ) %>%
  rename(height = bin_mid)

# --- 1.4 Plot the expected curve ---

speciesRichnessPlotUnbound <- ggplot(expected_curve_unbound, aes(y = height)) +
  # Potential MDE curve + ribbon
  geom_path(aes(x = mean_potential_richness, color = "Potential MDE"),
            size = 1, linetype = "dashed") +
  geom_ribbon(aes(xmin = mean_potential_richness - sd_potential_richness,
                  xmax = mean_potential_richness + sd_potential_richness,
                  fill = "Potential ±SD"), alpha = 0.2) +

  # Realised MDE curve + ribbon
  geom_path(aes(x = mean_realised_richness, color = "Realised MDE"),
            size = 1, linetype = "dashed") +
  geom_ribbon(aes(xmin = mean_realised_richness - sd_realised_richness,
                  xmax = mean_realised_richness + sd_realised_richness,
                  fill = "Realised ±SD"), alpha = 0.2) +

  scale_color_manual(
    name = "MDE Curve",
    values = c("Potential MDE" = "red", "Realised MDE" = "blue")
  ) +
  scale_fill_manual(
    name = "MDE ±SD",
    values = c("Potential ±SD" = "pink", "Realised ±SD" = "lightblue")
  ) +
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

DirectoryPlots <- file.path("../../figs/a5_plots_test/v8_real_niches_original_model_light_hum_temp")
filename <- file.path(DirectoryPlots,
                      "Replicate_1_SpeciesPool_1_unbound_MDnull_cappedPotential_vs_cappedRealised_199.pdf")

pdf(filename)
print(speciesRichnessPlotUnbound)
dev.off()

# --- Add surface area to the plot ---

# - Compute surface area per height

mh <- readRDS(file.path(mh_path, mh_file))
surface_area <- mh[,,, 1]
srface_area_per_height <- apply(surface_area, 3, sum) # in m2

# Add surface area to expected_curve_unbound
expected_curve_unbound <- expected_curve_unbound %>%
  mutate(surface_area = srface_area_per_height)

# - Plot with surface area

# Determine a scaling factor to align surface_area with richness on the same axis
# We'll use a simple linear transformation
range_richness <- range(c(expected_curve_unbound$mean_potential_richness,
                          expected_curve_unbound$mean_realised_richness))
range_surface <- range(expected_curve_unbound$surface_area)

# Linear scaling function to match surface_area to richness scale
scale_factor <- diff(range_richness) / diff(range_surface)
offset <- range_richness[1] - range_surface[1] * scale_factor

speciesRichnessPlotUnbound2 <- speciesRichnessPlotUnbound +
  # Surface area curve
  geom_path(aes(x = srface_area_per_height * scale_factor + offset, y = height, color = "Surface Area"),
            size = 1, linetype = "solid") +
  scale_color_manual(
    name = "MDE Curve / Surface Area",
    values = c("Potential MDE" = "red",
               "Realised MDE" = "blue",
               "Surface Area" = "darkolivegreen3")
  ) +
  scale_x_continuous(
    name = "Species Richness",
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Surface Area (m²)")
  )

# Save
filename <- file.path(DirectoryPlots,
                      "Replicate_1_SpeciesPool_1_unbound_MDnull_capped_Potential_vs_Realised_sa_sum_199.pdf")
pdf(filename)
print(speciesRichnessPlotUnbound2)
dev.off()