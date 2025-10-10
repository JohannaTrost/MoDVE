library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)

# #################################################################################################
#                               Why do some species shift upward?                                 #
# #################################################################################################


base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")
species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc")

# Get stats
species_distr_stats <- species_distr %>%
  group_by(Scenario, ForestID, SpeciesPool, SpeciesID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop")

# Sort scenarios and arrange data
species_distr_stats <- species_distr_stats %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

# Compute position shift
shift <- species_distr_stats %>%
  dplyr::select(ForestID, Scenario, SpeciesPool, SpeciesID, TimeStep, Year, Position) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = Position
  ) %>%
  mutate(
    diff = `CC` - `No CC`
  ) %>%
  filter(!is.na(diff))

summary_shift <- shift %>%
  group_by(SpeciesPool, SpeciesID) %>%
  arrange(Year, .by_group = TRUE) %>%  # ensure data is sorted by Year
  reframe(
    initial_diff = mean(diff[Year == 2000], na.rm = TRUE),
    diff = mean(tail(diff[!is.na(diff)], 20), na.rm = TRUE),
    last_year_alive = max(Year[!is.na(diff)]),
  ) %>%
  filter(last_year_alive >= 2050)

# Get species niches
niches <- NULL
for (sp in 1:10) {
  sp_niches <- read_csv(file.path(base_dir, "a2_1", paste0("SpeciesPool", sp, ".csv")),
                        show_col_types = FALSE) %>%
    dplyr::select("SpeciesID", "OptimumLight", "OptimumTemp", "OptimumHum",
           "MinLight", "MaxLight", "MinTemp", "MaxTemp", "MinHum", "MaxHum", "MaximumMass", "GrowthRate") %>%
    mutate(SpeciesPool = sp)

  if (is.null(niches)) {
    niches <- sp_niches
  } else {
    niches <- rbind(niches, sp_niches)
  }
}

# Compute ranges
niches$RangeTemp <- niches$MaxTemp - niches$MinTemp
niches$RangeHum <- niches$MaxHum - niches$MinHum
niches$RangeLight <- niches$MaxLight - niches$MinLight

# Combine with shift data
shift_niches <- summary_shift %>%
  left_join(niches, by = c("SpeciesPool", "SpeciesID"))

# Compute species distances
cols_to_use <- c("OptimumTemp", "OptimumHum", "MinTemp", "MaxTemp", "MinHum", "MaxHum",
                 "RangeTemp", "RangeHum", "OptimumLight", "MinLight", "MaxLight", "RangeLight")

shift_niches_dist <- shift_niches %>%
  group_by(SpeciesPool) %>%
  mutate(across(all_of(cols_to_use),
                ~ {
                    col_values <- .
                    sapply(seq_along(col_values), function(i) {
                      mean(abs(col_values[i] - col_values[-i]))
                    })
                  },
                .names = "MeanDist_{.col}")) %>%
  ungroup()

# - Compute niche overlap per species within group
# Function to compute pairwise overlap for one variable
range_overlap <- function(min1, max1, min2, max2) {
  overlap <- pmax(0, pmin(max1, max2) - pmax(min1, min2))
  union <- pmax(max1, max2) - pmin(min1, min2)
  overlap / union
}

compute_niche_distance <- function(df) {
  species_ids <- df$SpeciesID

  # All species pairs (excluding self)
  pairs <- expand.grid(SpeciesID1 = species_ids, SpeciesID2 = species_ids) %>%
    filter(SpeciesID1 != SpeciesID2)

  # Join data
  pairs <- pairs %>%
    left_join(df, by = c("SpeciesID1" = "SpeciesID")) %>%
    rename_with(~ paste0(., "_1"), -SpeciesID1) %>%
    left_join(df, by = c("SpeciesID2_1" = "SpeciesID")) %>%
    rename_with(~ paste0(., "_2"), -SpeciesID2_1)

  # Compute overlaps
  pairs <- pairs %>%
    mutate(
      overlap_temp  = range_overlap(MinTemp_1_2, MaxTemp_1_2, MinTemp_2, MaxTemp_2),
      overlap_hum   = range_overlap(MinHum_1_2, MaxHum_1_2, MinHum_2, MaxHum_2),
      overlap_light = range_overlap(MinLight_1_2, MaxLight_1_2, MinLight_2, MaxLight_2)
    )

  # Columns for which to compute absolute differences
  diff_cols <- c("OptimumTemp", "OptimumHum", "MinTemp", "MaxTemp", "MinHum", "MaxHum", "RangeTemp", "RangeHum",
                 "OptimumLight", "MinLight", "MaxLight", "RangeLight"
  )

  # Add average absolute differences dynamically
  for (col in diff_cols) {
    col1 <- paste0(col, "_1_2")
    col2 <- paste0(col, "_2")
    newcol <- paste0("AvgAbsDiff_", col)
    pairs[[newcol]] <- abs(pairs[[col1]] - pairs[[col2]])
  }

  # Summarise by focal species
  pairs <- pairs %>%
    group_by(SpeciesID1_2) %>%
    summarise(
      MaxOverlapTemp  = max(overlap_temp, na.rm = TRUE),
      MaxOverlapHum   = max(overlap_hum, na.rm = TRUE),
      MaxOverlapLight = max(overlap_light, na.rm = TRUE),
      across(starts_with("AvgAbsDiff_"), ~ mean(.x, na.rm = TRUE))
    ) %>%
    rename(SpeciesID = SpeciesID1_2)

  pairs
}

# Apply per ForestID + SpeciesPool
niche_overlap <- shift_niches %>%
  group_by(SpeciesPool) %>%
  group_modify(~ compute_niche_distance(.x)) %>%
  ungroup()

shift_niches_dist <- shift_niches_dist %>%
  left_join(., niche_overlap, by = c("SpeciesID", "SpeciesPool"))

# --- PCA

pca_data <- shift_niches_dist %>%
  select(MeanDist_OptimumHum, MaxOverlapHum,
         MeanDist_OptimumTemp, MaxOverlapTemp,
         MeanDist_OptimumLight, MaxOverlapLight,
         MaximumMass, GrowthRate
  ) %>%
  scale(.) %>%
  data.frame(.) %>%
  tibble(.)
res.pca <- dudi.pca(pca_data, scannf = FALSE, nf = 6)

# PCA results
DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/functional_analysis/MaxOverlap")
dir.create(DirectoryPlots)

# Variance explained
pdf(file.path(DirectoryPlots, "niche_dist_pca_scree.pdf"))
fviz_eig(res.pca)
dev.off()

pdf(file.path(DirectoryPlots, "niche_dist_pca_contrib_vars.pdf"))
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
dev.off()

pdf(file.path(DirectoryPlots, "niche_dist_pca_biplot_vars.pdf"))
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )
dev.off()

# -- Plot shift

pdf(file.path(DirectoryPlots, "niche_distance_pca_biplot_shift_traits.pdf"))
fviz_pca_biplot(res.pca,
                label = "var",
                labelsize = 3,
             #select.ind = list(ind = survived_idx),
             habillage = shift_niches_dist$diff > 0,
             palette = c("#00AFBB", "#FC4E07"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Shift",
             repel = TRUE)
dev.off()

# -- Plot denisty only

# Extract PCA coordinates for individuals
pca_coords <- as.data.frame(res.pca$li[, 1:2])  # First two PCs
colnames(pca_coords) <- c("PC1", "PC2")

# Add shift information
pca_coords <- pca_coords %>%
  mutate(diff = shift_niches_dist$diff) %>%
  mutate(Shift = ifelse(diff > 0, "Upward", "Downward"))

# Remove missing values (if any)
pca_coords <- na.omit(pca_coords)

# Compute group centroids
centroids <- pca_coords %>%
  group_by(Shift) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE)
  )

# Compute explained variance (percent)
eig_var <- res.pca$eig / sum(res.pca$eig) * 100

pdf(file.path(DirectoryPlots, "niche_dist_pca_shift_traits_group_centroids_density.pdf"), width = 7, height = 6)

ggplot(pca_coords, aes(x = PC1, y = PC2, color = Shift, fill = Shift)) +
  # Convex hulls
  geom_polygon(
    data = do.call(rbind, lapply(split(pca_coords, pca_coords$Shift), function(df) {
      df[chull(df$PC1, df$PC2), ]
    })),
    aes(group = Shift),
    alpha = 0.15,
    color = NA
  ) +
  # Density contours
  geom_density_2d(alpha = 0.3) +
  # Centroids
  geom_point(data = centroids, aes(x = PC1, y = PC2), shape = 21, size = 5, fill = "white", color = "black") +
  geom_text(data = centroids, aes(label = Shift), vjust = -1, fontface = "bold") +
  # Axis labels with variance explained
  labs(
    x = paste0("PC1 (", round(eig_var[1], 1), "%)"),
    y = paste0("PC2 (", round(eig_var[2], 1), "%)"),
    color = "Shift",
    fill = "Shift"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

dev.off()

# -- Plot mortality categories



# Below old

# --- Try MEM

pdf(file.path(DirectoryPlots, "diag_shift_hist.pdf"))
hist(shift_niches_dist$diff)
dev.off()

shapiro.test(shift_niches_dist$diff)

# Scale predictors
pca_data$diff <- shift_niches_dist$diff
lm_shift_range <- lm(
   diff ~ MeanDist_OptimumHum + MaxOverlapHum + MeanDist_OptimumTemp + MaxOverlapTemp +
    MeanDist_OptimumLight + MaxOverlapLight + MaximumMass + GrowthRate,
  data = pca_data)

summary(lm_shift_range)

# stepwise simplification
lm_shift_range_simpl <- lm(
   diff ~ MeanDist_OptimumHum + MaxOverlapTemp + GrowthRate,
  data = pca_data)

summary(lm_shift_range_simpl)


# - Try Binary response
pca_data$diff_bin <- as.numeric(shift_niches_dist$diff > 0)
glm_shift_bin <- glm(
  diff_bin ~ MeanDist_OptimumHum + MaxOverlapHum + MeanDist_OptimumTemp + MaxOverlapTemp +
    MeanDist_OptimumLight + MaxOverlapLight + MaximumMass + GrowthRate,
  data = pca_data,
  family = binomial(link = "logit")
)
summary(glm_shift_bin)

# -> Conclusion
# Species with higher growth rates are significantly less likely to shift upwards
# compared to those with lower growth rates (odds ratio = 0.30, p = 0.015).

# Relate to the number of recrutits - potential TODO

