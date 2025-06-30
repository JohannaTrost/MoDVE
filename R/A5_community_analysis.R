source("utils.R")

library(readr)
library(tidyr)
library(vegan)
library(dplyr)
library(purrr)
library(ggplot2)
library(forcats)
library(tidyverse)

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
  adjacentTurnover <- sapply(1:(nrow(betaMatrix) - 1), function(i) {
    betaMatrix[i, i + 1]
  })

  return(list(adjacentTurnover = adjacentTurnover,
              speciesHeightMatrix = speciesHeightMatrix))
}

DirectoryModelResults <- "tests/data/output_a4/"
DirectoryPlots <- "../../../figs/a5_plots_test/"
numSpeciesPools <- c(1, 2)
replicatePerSpeciesPool <- 1
timeStepStart <- 2
timeStepEnd <- 3

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
    for (timeStep in int_seq(timeStepStart, timeStepEnd)) {

      modelResultsPath <- paste0(DirectoryModelResults,
                                 "ID_SpeciesP_", numSpeciesPool, "_Rep_", rep,
                                 "/IndividualMatrixTimeStep", timeStep, ".csv")
      if (!file.exists(modelResultsPath)) {
        message("File does not exist: ", modelResultsPath)
        next
      }

      # Load results
      res <- read_csv(modelResultsPath)
      res <- res[res$Status == 1,]  # We are only interested in individuals that survived

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

      # Add to the turnover tibble
      turnoverDf <- turnoverDf %>%
              add_row(SpeciesPool = numSpeciesPool,
                      Replicate = rep,
                      TimeStep = timeStep,
                      LowerHeight = speciesHeightMatrix$height[-length(speciesHeightMatrix$height)],
                      UpperHeight = speciesHeightMatrix$height[-1],
                      BetaDivTurnover = adjacentTurnover,
                      Richness = )

    }

    # 1. Color coded plot of speciesHeightMatrix (height on x-axis, SpeciesID on y-axis, color = abundance)
    #    Last time step only
    pdf(paste0(DirectoryPlots, "SpeciesPool_", numSpeciesPool, "_Replicate_", rep, "_SpeciesAbundance.pdf"))
    print(PlotSpeciesHeightAbundance(res))
    dev.off()
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

# mkdir if not exists
if (!dir.exists(DirectoryPlots)) {
  dir.create(DirectoryPlots)
}

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
  pdf(filename, width = 6, height = 5)
  print(p)
  dev.off()
}

# --- 3. For each Replicate and SpeciesPool, make a line plot with height on y-axis and turnover on y-axis
#    each time step shouls have one line and it should have a color gradient (e.g. shapes of blue light to dark)

turnoverPlot <- turnoverDf %>%
  # Calculate midpoint height for plotting
  mutate(MidHeight = (LowerHeight + UpperHeight) / 2) %>%
  ggplot(aes(x = BetaDivTurnover, y = MidHeight, color = factor(TimeStep))) +
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

# --- 4. Richness per height plot

speciesRichnessPlot <- speciesHeight %>%
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