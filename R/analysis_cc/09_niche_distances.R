library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)
library(randomForest)
library(tibble)

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
  arrange(Year, .by_group = TRUE) %>%
  summarise(
    initial_diff = mean(.data$diff[.data$Year == 2000], na.rm = TRUE),
    valid_diff = .data$diff[!is.na(.data$diff)],
    diff = mean(tail(valid_diff, 20), na.rm = TRUE),
    diff_sem = sd(tail(valid_diff, 20), na.rm = TRUE) / sqrt(length(tail(valid_diff, 20))),
    last_year_alive = max(.data$Year[!is.na(.data$diff)]),
    .groups = "drop"
  ) %>%
  filter(last_year_alive >= 2050) %>%
  select(-valid_diff) %>%
  distinct()

write_csv(summary_shift, file.path(base_dir, "a5_species_shift_cc_vs_no_cc_last20_yrs_past_2050.csv"))

# Get species niches
niches <- NULL
for (sp in 1:10) {
  sp_niches <- read_csv(file.path(base_dir, "a2_1", paste0("SpeciesPool", sp, ".csv")),
                        show_col_types = FALSE) %>%
    dplyr::select("SpeciesID", "OptimumLight", "OptimumTemp", "OptimumHum",
                  "MinLight", "MaxLight", "MinTemp", "MaxTemp", "MinHum", "MaxHum", "MaximumMass",
                  "GrowthRate", "DispersalKernel"
    ) %>%
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

shift_species_q75 <- shift_niches_dist %>%
  filter(diff >= quantile(diff, 0.75))
write_csv(shift_species_q75,
          file.path(base_dir, "a5_species_shift_upwards_cc_vs_no_cc.csv"))

# --- PCA

pca_data <- shift_niches_dist %>%
  select(MeanDist_OptimumHum, MaxOverlapHum,
         MeanDist_OptimumTemp, MaxOverlapTemp,
         MeanDist_OptimumLight, MaxOverlapLight,
         MaximumMass, GrowthRate, MaxHum, MaxTemp, MinTemp
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



# --- Try MEM

pdf(file.path(DirectoryPlots, "diag_shift_hist.pdf"))
hist(shift_niches_dist$diff)
dev.off()

shapiro.test(shift_niches_dist$diff)

# Scale predictors
lm_data <- pca_data
lm_data$diff <- shift_niches_dist$diff
lm_data$OptimumTemp <- shift_niches_dist$OptimumTemp
lm_data$MinTemp <- shift_niches_dist$MinTemp
lm_data$MaxTemp <- shift_niches_dist$MaxTemp
lm_data$MaxHum <- shift_niches_dist$MaxHum
lm_data$DispersalKernel <- shift_niches_dist$DispersalKernel
lm_data$RangeTemp <- shift_niches_dist$RangeTemp
lm_shift_range <- lm(
   diff ~ MeanDist_OptimumHum + MaxOverlapHum + MeanDist_OptimumTemp + MaxOverlapTemp +
    MeanDist_OptimumLight + MaxOverlapLight + MaximumMass + GrowthRate + MaxTemp +
    MinTemp + MaxHum + DispersalKernel + OptimumTemp,
  data = lm_data)

summary(lm_shift_range)



# stepwise simplification
lm_shift_range_simpl <- lm(
   diff ~ MeanDist_OptimumHum * MaxOverlapTemp * GrowthRate,
  data = lm_data)

summary(lm_shift_range_simpl)

cor.test(lm_data$MaxOverlapTemp, lm_data$diff)
cor.test(lm_data$MaxTemp, lm_data$diff)
cor.test(lm_data$MaxTemp, lm_data$MaxOverlapTemp)

# - Try Binary response
lm_data$diff_bin <- as.factor(as.numeric(shift_niches_dist$diff > 0))
levels(lm_data$diff_bin)

glm_shift_bin <- glm(
  diff_bin ~ MeanDist_OptimumHum + MaxOverlapHum + MeanDist_OptimumTemp + MaxOverlapTemp +
    MeanDist_OptimumLight + MaxOverlapLight + MaximumMass + GrowthRate + MaxTemp,
  data = lm_data,
  family = binomial(link = "logit")
)
summary(glm_shift_bin)

# -> Conclusion
# Species with higher growth rates are significantly less likely to shift upwards
# compared to those with lower growth rates (odds ratio = 0.30, p = 0.0234).

glm_shift_bin <- glm(
  diff_bin ~ MaxOverlapHum + MaxOverlapTemp + GrowthRate + MaxTemp,
  data = lm_data,
  family = binomial(link = "logit")
)
summary(glm_shift_bin)

cor.test(lm_data$RangeTemp, lm_data$diff)

# - Make some plots

plt_data <- shift_niches_dist
plt_data$Shift <- ifelse(shift_niches_dist$diff > 0, "Upward", "Downward")

DirectoryPlots <- "../../figs/a5_plots_test/cc_vs_no_cc/position_shift"

# GROwth rate
pdf(file.path(DirectoryPlots, "Gr_boxplot_by_shift_final_species.pdf"), width = 10, height = 8)
ggplot(plt_data, aes(x = Shift, y = GrowthRate, fill = Shift)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Shift Category",
    y = "Growth rate",
    fill = "Species position shift"
  ) +
  theme_minimal()
dev.off()

# Max Temp
pdf(file.path(DirectoryPlots, "MaxTemp_boxplot_by_shift_final_species.pdf"), width = 10, height = 8)
ggplot(plt_data, aes(x = Shift, y = MaxTemp, fill = Shift)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Shift Category",
    y = "Maximum Temperature (°C)",
    fill = "Species position shift"
  ) +
  theme_minimal()
dev.off()

# Max temp overlap
pdf(file.path(DirectoryPlots, "MaxOverlapTemp_boxplot_by_shift_final_species.pdf"), width = 10, height = 8)
ggplot(plt_data, aes(x = Shift, y = MaxOverlapTemp, fill = Shift)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Shift Category",
    y = "Maximum temperature niche overlap (proportion)",
    fill = "Species position shift"
  ) +
  theme_minimal()
dev.off()

# Temperature range
pdf(file.path(DirectoryPlots, "RangeTemp_boxplot_by_shift_final_species.pdf"), width = 10, height = 8)
ggplot(plt_data, aes(x = Shift, y = RangeTemp, fill = Shift)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Shift Category",
    y = "Temperature range (°C)",
    fill = "Species position shift"
  ) +
  theme_minimal()
dev.off()

# ---- RF analysis

## set the seed to make your partition reproducible
set.seed(123)

down <- plt_data %>% filter(Shift == "Downward") %>% select(-diff, -last_year_alive, -SpeciesID, -GrowthRate)
up <- plt_data %>% filter(Shift == "Upward") %>% select(-diff, -last_year_alive, -SpeciesID, -GrowthRate)

n_down <- nrow(down)
n_up <- nrow(up)

train_ind_down <- sample(seq_len(n_down), size = floor(0.75 * n_down))
train_ind_up <- sample(seq_len(n_up), size = floor(0.75 * n_up))

train <- rbind(down[train_ind_down,], up[train_ind_up,])
test <- rbind(down[-train_ind_down,], up[-train_ind_up,])

sample_size <- min(table(train$Shift))

rf_model <- randomForest(
  factor(Shift) ~ ., data = train,
  importance = TRUE, ntree = 1000,
  sampsize = c('Downward'= sample_size,'Upward'= sample_size)
)

# --- 3. Model performance ----------------------------------------------

print(rf_model)
# Look at out-of-bag error rate and confusion matrix

prediction <- data.frame(predict(rf_model, test, type='prob'))
preds <- ifelse(prediction$Downward > prediction$Upward, "Downward", "Upward")
# ensure test$Shift has the same levels in the same order
test$Shift <- factor(test$Shift, levels = c("Downward", "Upward"))
# make preds an ordered factor with the same levels
preds_ord <- factor(preds, levels = c("Downward", "Upward"), ordered = TRUE)

# compute ROC
roc1 <- roc(test$Shift, preds_ord)
print(roc1)

# --- 4. Variable importance --------------------------------------------

importance_df <- as.data.frame(importance(rf_model))
importance_df <- importance_df %>%
  rownames_to_column("Trait") %>%
  arrange(desc(MeanDecreaseGini))

# Print ranked trait importance
print(importance_df)

# --- 5. Visualize variable importance ----------------------------------

pdf(file.path(DirectoryPlots, "randomForestΤraitImportanceShift.pdf"), width = 7, height = 7)
ggplot(importance_df, aes(x = reorder(Trait, MeanDecreaseGini),
                          y = MeanDecreaseGini)) +
  geom_col(fill = "#56B4E9") +
  coord_flip() +
  labs(x = "Trait",
       y = "Mean Decrease in Gini (importance)") +
  theme_minimal(base_size = 13)
dev.off()

# --- Partial dependence (to interpret top traits) ----------

# For the most important traits:
library(pdp)

top_traits <- importance_df$Trait[1:3]
pred_df <- as.data.frame(test)

for (trait in top_traits) {
  outfile <- file.path(DirectoryPlots,
                       paste0("randomForest_", trait, "_PartialDepShiftUpTrain.pdf"))
  pdf(outfile, width = 7, height = 5)
    # Get partial dependence data
  pd <- partial(
    object = rf_model,
    pred.var = trait,
    which.class = "Upward",
    prob = TRUE,
    train = pred_df
  )

  # Create the plot
  print(autoplot(pd))
  dev.off()
}


# ---- REDUCED model removing correlated predictors

# Load required package
library(caret)

# --- Step 1: Define response and predictors ---
response <- "Shift"
predictors <- setdiff(names(plt_data), response)

# --- Step 2: Subset predictors only ---
predictor_data <- plt_data[, predictors]

# --- Step 3: Calculate correlation matrix ---
cor_matrix <- cor(predictor_data, use = "pairwise.complete.obs")

# --- Step 4: Identify highly correlated predictors ---
# Adjust cutoff (e.g., 0.7, 0.8, or 0.9 depending on strictness)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE)

# --- Step 5: Remove correlated predictors ---
plt_data_reduced <- plt_data[, !(names( ) %in% high_cor)]

# --- Step 6: Check remaining predictors ---
cat("Removed predictors due to high correlation:\n")
print(high_cor)
