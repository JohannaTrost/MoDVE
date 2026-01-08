source("utils.R")

library(readr)
library(tidyr)
library(vegan)
library(dplyr)
library(purrr)
library(ggplot2)
library(forcats)
library(tidyverse)
library(patchwork)
library(zoo)

PlotSpeciesHeightAbundance <- function(res) {
  abundance_df <- res %>%
  group_by(height, SpeciesID) %>%
  summarise(abundance = n(), .groups = "drop")

  # Convert SpeciesID to factor for proper x-axis ordering
  abundance_df <- abundance_df %>%
    mutate(SpeciesID = factor(SpeciesID, levels = sort(unique(SpeciesID))))

  plt <-ggplot(abundance_df, aes(x = SpeciesID, y = height, fill = abundance)) +
    geom_tile(color = "grey70") +
    scale_fill_gradient(low = "grey", high = "black") +
    theme_minimal() +
    labs(x = "Species", y = "Height", fill = "Abundance") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid = element_blank())

  return(plt)
}

ComputeDivMetrics <- function(data) {
  richness <- length(unique(data$SpeciesID))
  abundance <- dim(data)[1]
  return(
      list(
        richness = richness,
        abundance = abundance,
        mass = sum(data$Mass),
        evenness = ComputeEvenness(data, abundance, richness)
      )
  )
}

ComputeEvenness <- function(data, totalAbundance, richness) {
  # Count rows per SpeciesID
  speciesStats <- data %>%
    group_by(SpeciesID) %>%
    summarise(Abundance = n()) %>%
    ungroup() %>%
    mutate(Proportion = Abundance / totalAbundance)

  simpsons_d <-  1 / sum(speciesStats$Proportion^2) # Probability of two individuals being from the same species
  evenness <- simpsons_d / richness # Simpson's evenness

  return(evenness)
}

ComputeDivTurnover <- function(data) {
  # 1. Compute Jaccard dissimilarity
  speciesHeightMatrix <- data %>%
    group_by(height, SpeciesID) %>%
    summarise(present = 1, .groups = "drop") %>%
    pivot_wider(names_from = SpeciesID, values_from = present, values_fill = 0) %>%
    arrange(height)
  betaMatrix <- speciesHeightMatrix %>%
    select(-height) %>%
    vegdist(., method = "jaccard", binary = TRUE) %>%
    as.matrix()

  # 2. Get turnover between adjacent heights
  if (nrow(betaMatrix) > 1) {
    adjacentTurnover <- sapply(1:(nrow(betaMatrix) - 1), function(i) {
      betaMatrix[i, i + 1]
    })
  } else {
    adjacentTurnover <- numeric(0)  # No adjacent heights to compare
  }


  return(list(adjacentTurnover = adjacentTurnover,
              speciesHeightMatrix = speciesHeightMatrix))
}

DirectoryModelResults <- "/Users/johanna/Uni/masterarbeit/data/modve_output/regua/climdata_era5_cmip6_1981-2100_ssp245/a5_upward_shifted/forest0/"
DirectoryPlots <- "../../../figs/a5_plots_test/climdata_era5_cmip6_1981-2100_ssp245/forest0/a5_upward_shifted/"
numSpeciesPools <- seq(1, 10)
replicatePerSpeciesPool <- 1
timeStepStart <- 100
timeStepEnd <- 199
stepSize <- 5

# mkdir if not exists
if (!dir.exists(DirectoryPlots)) {
  dir.create(DirectoryPlots, recursive = TRUE)
}

divMetricsDf <- tibble(
  SpeciesPool = integer(),
  Replicate = integer(),
  TimeStep = integer(),
  Richness = numeric(),
  Abundance = numeric(),
  Mass = numeric(),
  Evenness = numeric()
)
turnoverDf <- tibble(
  SpeciesPool = integer(),
  Replicate = integer(),
  TimeStep = integer(),
  LowerHeight = numeric(),
  UpperHeight = numeric(),
  BetaDivTurnover = numeric(),
)
speciesHeight <- tibble(
  SpeciesPool = integer(),
  Replicate = integer(),
  TimeStep = integer(),
  SpeciesID = integer(),
  height = numeric(),
)

for (rep in seq_len(replicatePerSpeciesPool)) {
  for (numSpeciesPool in numSpeciesPools) {

    # --- 5. Analysis of mortality and recruits
    ComSum <- read_csv(paste0(DirectoryModelResults, "ID_SpeciesP_", numSpeciesPool, "_Rep_", rep, "/CommunitySummary.csv"))

    # Reshape data to long format
    ComSumLong <- ComSum %>%
      select(timeStep, Recruits, starts_with("Mortality")) %>%
      pivot_longer(
        cols = -timeStep,
        names_to = "Type",
        values_to = "Count"
      )

    # Create the ggplot
    DemogPlot <- ggplot(ComSumLong, aes(x = timeStep, y = Count, color = Type)) +
      geom_line(alpha = 0.5) +
      labs(
        x = "Time Step",
        y = "Count",
        color = "Event Type"
      ) +
      theme_minimal()

    filename <- paste0(DirectoryPlots, "SpeciesPool_", numSpeciesPool, "_Replicate_", rep, "_Demography.pdf")
    pdf(filename, width = 8, height = 4)
    print(DemogPlot)
    dev.off()

    # --- Collecting data for each time step
    for (timeStep in int_seq(timeStepStart, timeStepEnd)) {

      modelResultsPath <- paste0(DirectoryModelResults,
                                 "ID_SpeciesP_", numSpeciesPool, "_Rep_", rep,
                                 "/IndividualMatrixTimeStep", timeStep, ".csv")
      if (!file.exists(modelResultsPath)) {
        message("File does not exist: ", modelResultsPath)
        next
      }
      print(timeStep)
      # Load results
      res <- read_csv(modelResultsPath)
      res <- res[res$Status == 1,]  # We are only interested in individuals that survived

      if (nrow(res) == 0) {
        message("No individuals survived at time step ", timeStep, " for Species Pool ", numSpeciesPool, " and Replicate ", rep)
        next
      }

      # Collect species height data
      res["height"] <- res$Z - 0.5
      speciesHeight <- speciesHeight %>%
        add_row(SpeciesPool = numSpeciesPool,
                Replicate = rep,
                TimeStep = timeStep,
                SpeciesID = res$SpeciesID,
                height = res$height)

      # Compute different metrics to assess the community diversity and productivity
      speciesDiv <- ComputeDivMetrics(res)

      # Add to the dataframe
      divMetricsDf <- divMetricsDf %>%
      add_row(SpeciesPool = numSpeciesPool,
              Replicate = rep,
              TimeStep = timeStep,
              Richness = speciesDiv$richness,
              Abundance = speciesDiv$abundance,
              Mass = speciesDiv$mass,
              Evenness = speciesDiv$evenness)

      # - Diversity turnover over heights
      divTurnover <- ComputeDivTurnover(res)
      speciesHeightMatrix <- divTurnover$speciesHeightMatrix
      adjacentTurnover <- divTurnover$adjacentTurnover

      print(timeStep)

      # Add to the turnover tibble
      turnoverDf <- turnoverDf %>%
              add_row(SpeciesPool = numSpeciesPool,
                      Replicate = rep,
                      TimeStep = timeStep,
                      LowerHeight = speciesHeightMatrix$height[-length(speciesHeightMatrix$height)],
                      UpperHeight = speciesHeightMatrix$height[-1],
                      BetaDivTurnover = adjacentTurnover)

    }

    if (nrow(res) > 0) {
      # 1. Color coded plot of speciesHeightMatrix (height on x-axis, SpeciesID on y-axis, color = abundance)
      #    Last time step only
      abundancePlotPath <- paste0(DirectoryPlots, "SpeciesPool_", numSpeciesPool, "_Replicate_", rep, "_SpeciesAbundance.pdf")
      pdf(abundancePlotPath)
      print(PlotSpeciesHeightAbundance(res))
      dev.off()
      print(abundancePlotPath)
    }
  }
}

# Species height distributions
heightDistr <- speciesHeight %>% group_by(Replicate, SpeciesPool, SpeciesID) %>%
  summarise(meanHeight = mean(height, na.rm = TRUE),
            p05 = quantile(height, 0.05, na.rm = TRUE),
            p25 = quantile(height, 0.25, na.rm = TRUE),
            p50 = quantile(height, 0.50, na.rm = TRUE),
            p75 = quantile(height, 0.75, na.rm = TRUE),
            p95 = quantile(height, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

# ------ Visualization ------

# --- 2. Box plots x = SpeciesID y = height for each SpeciesPool and Replicate

# Split data by SpeciesPool and Replicate
split_data <- speciesHeight %>%
  mutate(
    SpeciesID = as.factor(SpeciesID),
    SpeciesPool = as.factor(SpeciesPool),
    Replicate = as.factor(Replicate)
  ) %>%
  group_split(SpeciesPool, Replicate)

# Assuming DirectoryPlots is a string with a trailing slash, e.g. "outputs/plots/"
for (df in split_data) {
  # Extract identifiers
  sp <- unique(df$SpeciesPool)
  rep <- unique(df$Replicate)

  # Reorder species by median height
  df <- df %>%
    mutate(SpeciesID = fct_reorder(SpeciesID, height, .fun = median))

  # Create the plot
  p <- ggplot(df, aes(x = SpeciesID, y = height, fill = SpeciesPool)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black") +
    labs(
      x = "Species ID",
      y = "Height (m)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Set3")

  # File name
  filename <- paste0(DirectoryPlots, "SpeciesPool_", sp, "_Replicate_", rep, "_HeightDistribution.pdf")

  # Save plot to PDF
  pdf(filename, width = 18, height = 9)
  print(p)
  dev.off()
}

# --- 3. For each Replicate and SpeciesPool, make a line plot with height on y-axis and turnover on y-axis
#    each time step shouls have one line and it should have a color gradient (e.g. shapes of blue light to dark)

turnoverPlot <- turnoverDf %>%
  filter(TimeStep %% 5 == 0) %>%
  # Calculate midpoint height for plotting
  mutate(MidHeight = (LowerHeight + UpperHeight) / 2) %>%
  ggplot(aes(x = BetaDivTurnover, y = MidHeight, color = factor(TimeStep))) +
  geom_path(aes(group = TimeStep), linewidth = 1) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    name = "Time Step",
    values = colorRampPalette(c("lightblue", "darkblue"))(length(unique(turnoverDf$TimeStep)))
  ) +
  # Create facets for Replicate x SpeciesPool
  facet_grid(Replicate ~ SpeciesPool,
             labeller = labeller(
               Replicate = function(x) paste("Replicate", x),
               SpeciesPool = function(x) paste("Species Pool", x)
             )) +
  # Customize labels and theme
  labs(
    x = "Beta Diversity Turnover",
    y = "Height (m)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

filename <- paste0(DirectoryPlots, "SpeciesPool_", sp, "_Replicate_", rep, "_VerticalDivTurnover.pdf")

# Save plot to PDF
pdf(filename)
print(turnoverPlot)
dev.off()

# -- Turnover of final time steps (avg/sd)

ts <- 160

agg_data <- turnoverDf %>%
  filter(TimeStep >= ts) %>%
  group_by(SpeciesPool, LowerHeight, UpperHeight) %>%
  summarise(
    BetaDivTurnover_mean = mean(BetaDivTurnover, na.rm = TRUE),
    BetaDivTurnover_low = quantile(BetaDivTurnover, 0.05, na.rm = TRUE),
    BetaDivTurnover_high = quantile(BetaDivTurnover, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(MidHeight = (LowerHeight + UpperHeight) / 2)

turnoverPlot_final <- agg_data %>%
  ggplot(aes(y = MidHeight)) +
  geom_ribbon(
    aes(
      xmin = BetaDivTurnover_low,
      xmax = BetaDivTurnover_high
    ),
    fill = "lightblue",
    alpha = 0.3
  ) +
  geom_path(
    aes(x = BetaDivTurnover_mean),
    color = "darkblue",
    linewidth = 1
  ) +
  geom_point(
    aes(x = BetaDivTurnover_mean),
    size = 1, alpha = 0.8, color = "darkblue"
  ) +
  facet_grid(
    . ~ SpeciesPool,
    labeller = labeller(
      SpeciesPool = function(x) paste("Species Pool", x)
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      ceiling(max(agg_data$MidHeight)),
      by = 5
    )
  ) +
  labs(
    x = "Beta Diversity Turnover",
    y = "Height (m)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 12),       # tick label size
    axis.title = element_text(size = 14),      # axis titles larger
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

filename <- paste0(DirectoryPlots, "VerticalDivTurnover_", ts, ".pdf")
pdf(filename)
print(turnoverPlot_final)
dev.off()

# --- 4. Richness per height plot

speciesRichnessPlot <- speciesHeight %>%
  filter(TimeStep %% 5 == 0) %>%
  group_by(Replicate, SpeciesPool, TimeStep, height) %>%
  summarise(Richness = n_distinct(SpeciesID), .groups = "drop") %>%
  ggplot(aes(x = Richness, y = height, color = factor(TimeStep))) +
  geom_path(aes(group = TimeStep), size = 1) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    name = "Time Step",
    values = colorRampPalette(c("lightblue", "darkblue"))(length(unique(turnoverDf$TimeStep)))
  ) +
  # Create facets for Replicate x SpeciesPool
  facet_grid(Replicate ~ SpeciesPool,
             labeller = labeller(
               Replicate = function(x) paste("Replicate", x),
               SpeciesPool = function(x) paste("Species Pool", x)
             )) +
  # Customize labels and theme
  labs(
    x = "Species Richness",
    y = "Height (m)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

filename <- paste0(DirectoryPlots, "SpeciesPool_", sp, "_Replicate_", rep, "_VerticalSpeciesRichness.pdf")

pdf(filename)
print(speciesRichnessPlot)
dev.off()

# Final richness plot -> Showing last time step
ts <- 199


# --- 4. Richness and abundance over time (overall and per height bins)

overallStats <- speciesHeight %>%
  group_by(TimeStep, Replicate, SpeciesPool) %>%
  summarize(
    Richness = length(unique(SpeciesID)),
    Abundance = n()
  )

# Assuming your tibble is named `df`
# Example: df <- tibble::tibble(TimeStep = ..., Richness = ..., Abundance = ...)

# Plot for Richness over Time
p1 <- ggplot(overallStats, aes(x = TimeStep, y = Richness,
                               group = interaction(SpeciesPool, Replicate),
                               color = as.factor(SpeciesPool))) +
  geom_line(size = 1) +
  labs(x = "Time Step", y = "Richness", color = "Species Pool") +
  theme_minimal()

# Plot for Abundance over Time
p2 <- ggplot(overallStats, aes(x = TimeStep, y = Abundance,
                               group = interaction(SpeciesPool, Replicate),
                               color = as.factor(SpeciesPool))) +
  geom_line(size = 1) +
  labs(x = "Time Step", y = "Abundance", color = "Species Pool") +
  theme_minimal()

# Combine plots with shared legend
combined_plot <- (p1 + p2) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# Save to PDF
filename <- paste0(DirectoryPlots, "OverallRichnessAbundance.pdf")
ggsave(filename, plot = combined_plot, width = 8, height = 4)

# --- Mid domain randomization to compare to actual richness

# --- Step 1. Extract realized vertical ranges ---
vertical_ranges <- speciesHeight %>%
  filter(TimeStep == 199 & SpeciesPool == 1 & Replicate == 1) %>%
  group_by(SpeciesID) %>%
  summarise(height_range = max(height) - min(height), .groups = "drop")

# Maximum observed height
max_height <- speciesHeight %>%
  filter(TimeStep == 199, SpeciesPool == 4, Replicate == 1) %>%
  summarise(max_height = max(height)) %>%
  pull(max_height)

bins <- 0:ceiling(max_height)

# --- Step 2. Function to run one MDR ---
run_mdr <- function() {
  vr <- vertical_ranges %>%
    mutate(max_height = max_height - height_range,
           min_h = runif(n(), min = 0, max = max_height),
           max_h = min_h + height_range)

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

# --- Step 3. Run many randomizations (Monte Carlo) ---
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

# --- Smooth actual richness data ---
smooth_data <- speciesHeight %>%
  filter(TimeStep == ts) %>%
  group_by(Replicate, SpeciesPool, TimeStep, height) %>%
  summarise(Richness = n_distinct(SpeciesID), .groups = "drop") %>%
  arrange(Replicate, SpeciesPool, height) %>%
  group_by(Replicate, SpeciesPool) %>%
  mutate(
    Richness_smoothed = rollapply(
      Richness,
      width = 7,
      FUN = mean,
      align = "center",
      fill = NA,
      partial = TRUE
    )
  )

# --- Combine actual and expected data for plotting ---
combined_sim_null <- smooth_data %>%
  filter(SpeciesPool == 1, Replicate == 1) %>%
  left_join(expected_curve, by = "height")

# --- Step 5. Plot the smoothed actual richness against expected curve ---
speciesRichnessPlot <- ggplot(combined_sim_null, aes(x = Richness_smoothed, y = height)) +
  geom_path(aes(color = "Actual Richness"), size = 1) +
  geom_point(aes(color = "Actual Richness"), size = 2, alpha = 0.8) +
  geom_path(aes(x = mean_richness, y = height, color = "Avg. MDE Richness"), size = 1, linetype = "dashed") +
  geom_ribbon(aes(y = height, xmin = mean_richness - sd_richness, xmax = mean_richness + sd_richness, fill = "MDE +-SD"), alpha = 0.2) +
  scale_color_manual(name = "Richness Type", values = c("Actual Richness" = "darkblue", "Avg. MDE Richness" = "red")) +
  scale_fill_manual(name = "MDE +-SD", values = c("MDE +-SD" = "lightgrey")) +
  facet_grid(Replicate ~ SpeciesPool,
             labeller = labeller(
               Replicate = function(x) paste("Replicate", x),
               SpeciesPool = function(x) paste("Species Pool", x)
             )) +
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

filename <- paste0(DirectoryPlots, "Replicate_1_SpeciesPool_1_MDnull_VerticalSpeciesRichness_", ts, ".pdf")

pdf(filename)
print(speciesRichnessPlot)
dev.off()

# --- Plot the smoothed actual richness against height ---
speciesRichnessPlot_final <- smooth_data %>%
  ggplot(aes(x = Richness_smoothed, y = height)) +
  #geom_path(aes(), size = 1, color = "darkblue",) +
  geom_path(aes(x = Richness_smoothed, y = height),
            color = "darkblue", linewidth = 1) +
  geom_point(size = 1, alpha = 0.8, color = "darkblue") +
  # Create facets for Replicate x SpeciesPool
  facet_grid(Replicate ~ SpeciesPool,
             labeller = labeller(
               Replicate = function(x) paste("Replicate", x),
               SpeciesPool = function(x) paste("Species Pool", x)
             )) +
  scale_y_continuous(
    breaks = seq(
      0,
      ceiling(max(smooth_data$height)),
      by = 5
    )
  ) +
  # Customize labels and theme
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

filename <- paste0(DirectoryPlots,"Replicate_", rep, "_SmoothVerticalSpeciesRichness_", ts, ".pdf")

pdf(filename)
print(speciesRichnessPlot_final)
dev.off()