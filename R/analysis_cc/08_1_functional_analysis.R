library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)
library(vegan)
library(tibble)
library(randomForest)
library(patchwork)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

# Get species niches
niches <- NULL
for (sp in 1:10) {
  sp_niches <- read_csv(file.path(base_dir, "a2_1", paste0("SpeciesPool", sp, ".csv")),
                        show_col_types = FALSE) %>%
    dplyr::select(-LightResponseA, -LightResponseB, -LightResponseC, -MinWind, -MaxWind, -OptimumWind,
           -DispersalKernelWindEffect) %>%
    mutate(SpeciesPool = sp)

  if (is.null(niches)) {
    niches <- sp_niches
  } else {
    niches <- rbind(niches, sp_niches)
  }
}

# Scale without SpeciesID and SpeciesPool and add back
niches_scaled <- niches %>%
  dplyr::select(-SpeciesID, -SpeciesPool) %>%
  mutate(across(everything(), ~ scale(.)[,1])) %>%
  bind_cols(niches %>% dplyr::select(SpeciesID, SpeciesPool), .)

# compute sds (NA if column is all NA)
sds <- sapply(niches_scaled, sd, na.rm = TRUE)
# columns with zero sd or NA sd
problem_cols <- names(sds)[is.na(sds) | sds == 0]

# PCA
pca_data <- niches_scaled %>%
  dplyr::select(-all_of(problem_cols)) %>%
  # - RecruitmentInvestmentRel, AgeAtMaturity
  dplyr::select(-MaxRecruitsAtMaxMass, - MassAtMaturity) %>%  # Remove cols with similar contribution
  mutate(
    RangeTemp = MaxTemp - MinTemp,
    RangeHum = MaxHum - MinHum,
    RangeLight = MaxLight - MinLight
  )

res.pca <- pca_data %>%
  dplyr::select(-SpeciesID, -SpeciesPool) %>%
  dudi.pca(., scannf = FALSE, nf = 10)

# PCA results
DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/functional_analysis")
dir.create(DirectoryPlots)

# Variance explained
pdf(file.path(DirectoryPlots, "trait_pca_scree.pdf"))
fviz_eig(res.pca)
dev.off()

pdf(file.path(DirectoryPlots, "trait_pca_contrib_vars.pdf"))
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
dev.off()

pdf(file.path(DirectoryPlots, "trait_pca_biplot_vars_v2.pdf"))
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )
dev.off()

# --- Species survival categories

species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

# Get stats
species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario) %>%
  dplyr::select(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year, Position) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = Position
  ) %>%
  mutate(
    diff = `CC` - `No CC`
  )

# Survival
species_survival_cat <- species_distr_stats %>%
  filter(Year >= 2080) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    has_NoCC = any(!is.na(`No CC`)),
    has_CC   = any(!is.na(CC)),
    .groups = "drop"
  ) %>%
  mutate(
    survival = case_when(
      has_CC & !has_NoCC ~ "Survived CC",
      !has_CC & has_NoCC ~ "Died with CC",
      has_CC & has_NoCC  ~ "Survived CC",
      TRUE ~ "none"  # optional: if both are NA
    )
  ) %>%
  dplyr::select(SpeciesPool, SpeciesID, survival)

# Merge with niches
pca_data_surv <- pca_data %>%
  left_join(., species_survival_cat, by = c("SpeciesPool", "SpeciesID")) %>%
  mutate(survival = ifelse(is.na(survival), "none", survival))

# Biplot with survival categories

# Add a grouping column with NAs for excluded individuals
survived_idx <- which(pca_data_surv$survival != "none")
pca_data_surv$survival_filtered <- pca_data_surv$survival
pca_data_surv$survival_filtered[-survived_idx] <- NA

pdf(file.path(DirectoryPlots, "trait_pca_biplot_survival_v2.pdf"))
fviz_pca_biplot(res.pca,
                label = "var",
                labelsize = 3,
             select.ind = list(ind = survived_idx),
             habillage = pca_data_surv$survival_filtered,
             palette = c("#00AFBB", "#FC4E07", "#EFD2CB", "#241623"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Survival Groups",
             repel = TRUE)
dev.off()

# --- Boxplots of traits by survival category

vars <- c("GrowthRate", "OptimumHum", "RangeTemp")

plt_data <- left_join(species_survival_cat, niches, by = c("SpeciesPool", "SpeciesID"))
plt_data$RangeTemp <- plt_data$MaxTemp - plt_data$MinTemp

for (var in vars) {
  pdf(file.path(DirectoryPlots, paste0("trait_boxplot_", var, "_survival.pdf")), width = 6, height = 4)
  print(
    ggplot(plt_data, aes_string(x = "survival", y = var, fill = "survival")) +
      geom_boxplot(notch = TRUE) +
      scale_fill_manual(values = c("#00AFBB", "#FC4E07", "#EFD2CB", "#241623")) +
      labs(x = "Survival Category", y = var) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
  )
  dev.off()
}

# --- 1. RF model data ---------------------------------

rf_data <- plt_data
rf_data$RangeHum <- rf_data$MaxHum - rf_data$MinHum
rf_data$RangeLight <- rf_data$MaxLight - rf_data$MinLight
rf_data$survival <- as.factor(rf_data$survival)
trait_cols <- names(plt_data)[!(names(plt_data) %in%
  c("SpeciesID", "survival"))]

# --- 2. RF model to predict survival --------------------------------

set.seed(123)

train_idx <- sample(seq_len(nrow(rf_data)), size = floor(0.75 * nrow(rf_data)))

train <- rf_data[train_idx,]
test <- rf_data[-train_idx,]

rf_model_surv <- randomForest(
  survival ~ ., data = rf_data[, c("survival", trait_cols)],
  importance = TRUE, ntree = 1000,
  sampsize = rep(min(table(rf_data$survival)), 2)
)

# --- 3. Model performance ----------------------------------------------

print(rf_model_surv)
# Look at out-of-bag error rate and confusion matrix -> 24.09%

prediction <- data.frame(predict(rf_model_surv, test, type='prob'))
preds <- ifelse(prediction$Died.with.CC > prediction$Survived.CC, "Died with CC", "Survived CC")
# ensure test$Shift has the same levels in the same order
test$survival <- factor(test$survival, levels = c("Died with CC", "Survived CC"))
# make preds an ordered factor with the same levels
preds_ord <- factor(preds, levels = c("Died with CC", "Survived CC"), ordered = TRUE)

# compute ROC
test <- test %>%
  mutate(
    survival_num = case_when(
      survival == "Died with CC" ~ 0,
      survival == "Survived CC"  ~ 1,
      TRUE ~ NA_real_
    )
  )
preds_ord_num <- recode(preds_ord,
  "Died with CC" = 0,
  "Survived CC" = 1
)
roc1 <- roc(test$survival_num, preds_ord_num, auc = TRUE)
print(roc1)

# --- 4. Variable importance --------------------------------------------

importance_df_surv <- as.data.frame(importance(rf_model_surv))
importance_df_surv <- importance_df_surv %>%
  rownames_to_column("Trait") %>%
  arrange(desc(MeanDecreaseGini))

# Print ranked trait importance
print(importance_df_surv)

# --- 5. Visualize variable importance below

# --- 6. Optional: partial dependence (to interpret top traits) ----------

# For the most important traits:
library(pdp)

top_traits <- importance_df$Trait[1:2]
pred_df <- as.data.frame(rf_data)

for (trait in top_traits) {
  outfile <- file.path(DirectoryPlots,
                       paste0("randomForest_", trait, "_PartialDepSurvivalnoCC.pdf"))
  pdf(outfile, width = 7, height = 5)
    # Get partial dependence data
  pd <- partial(
    object = rf_model_balanced,
    pred.var = trait,
    which.class = "only No CC",
    prob = TRUE,
    train = pred_df
  )

  # Create the plot
  print(autoplot(pd))
  dev.off()
}

# --------- Position shift

year <- 2080
species_shift <- species_distr_stats %>%
  filter(Year >= year) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    AvgDiff = mean(diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Shift = ifelse(AvgDiff > 0, "Upward", "Downward"))

pca_data_shift <- pca_data %>%
  left_join(., species_shift, by = c("SpeciesPool", "SpeciesID")) #%>%
  #mutate(UpwardShift = ifelse(is.na(UpwardShift), "none", mortality))

# Exclude non existant species
na_idx <- which(is.na(pca_data_shift$AvgDiff))

pdf(file.path(DirectoryPlots, paste0("trait_pca_biplot_shift_all_traits_", year, "-2100.pdf")))
fviz_pca_biplot(res.pca,
                label = "var",
                labelsize = 3,
             select.ind = list(ind = na_idx),
             habillage = pca_data_shift$Shift,
             palette = c("#00AFBB", "#FC4E07"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Vertical shift with CC",
             repel = TRUE)
dev.off()

# -- Plot denisty only

# Extract PCA coordinates for individuals
pca_coords <- as.data.frame(res.pca$li[, 1:2])  # First two PCs
colnames(pca_coords) <- c("PC1", "PC2")

# Add survival information
pca_coords$Shift <- pca_data_shift$Shift

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

pdf(file.path(DirectoryPlots, "trait_pca_shift_group_centroids_density.pdf"), width = 7, height = 6)

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


# Boxplot position by survival category
species_survival_cat <- species_survival_cat %>%
  mutate(survival = ifelse(survival == "both" | survival == "only CC", "Survived CC", "Survived without CC"))
species_survival_cat <- species_survival_cat %>%rename(Scenario = survival) %>%
  mutate(Scenario = ifelse(Scenario == "Survived CC", "CC", "No CC"))

species_pos <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario) %>%
  filter(year >= 2080) %>%
  group_by(SpeciesID, Scenario, SpeciesPool) %>%
  summarise(AvgPosition = mean(Position, na.rm = TRUE),
            .groups = "drop")

species_pos_surv <- inner_join(species_pos, species_survival_cat,
                               by = c("SpeciesPool", "SpeciesID", "Scenario")) %>%
  left_join(., niches, by = c("SpeciesID", "SpeciesPool")) %>%
  mutate(Scenario = ifelse(Scenario == "CC", "Survived CC", "Died with CC")) %>%
  mutate(Scenario = factor(Scenario))

# PLot survival and relationship between growth rate and position
pdf(file.path(DirectoryPlots, "gr_vs_pos_col_scenario.pdf"), width = 10, height = 8)
ggplot(species_pos_surv, aes(x = GrowthRate, y = AvgPosition, color = Scenario)) +
  geom_point() +
  #scale_color_viridis_d(option = "plasma") +
  labs(
    x = "Growth Rate",
    y = "Position (m)",
    color = "Survival"
  ) +
  theme_minimal()
dev.off()

# Do boxplot of AvgPosition by Scenario
pdf(file.path(DirectoryPlots, "boxplot_position_by_survival_CC.pdf"), width = 6, height = 4)
ggplot(species_pos_surv, aes(x = Scenario, y = AvgPosition, fill = Scenario)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Survival",
    y = "Average Position (m)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
dev.off()

wilcox.test(species_pos_surv[species_pos_surv$Scenario == "Survived CC", ]$GrowthRate,
            species_pos_surv[species_pos_surv$Scenario == "Died with CC", ]$GrowthRate)

pdf(file.path(DirectoryPlots, "boxplot_gr_by_survival_CC.pdf"), width = 6, height = 4)
ggplot(species_pos_surv, aes(x = Scenario, y = GrowthRate, fill = Scenario)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Survival",
    y = "GrowthRate"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
dev.off()

# ---- Random Forest for shift

# - look at correlations of variables

# Compute correlation matrix (exclude non-numeric columns automatically)
cor_matrix <- cor(pca_data_shift %>% dplyr::select(where(is.numeric)), use = "complete.obs")

# Extract correlations with AgeAtMaturity
cor_age <- cor_matrix[,"AgeAtMaturity"]

# Filter variables with high absolute correlation (> 0.7)
high_cor_vars <- names(cor_age[abs(cor_age) > 0.7 & names(cor_age) != "AgeAtMaturity"])

# Print results
cat("Variables highly correlated with AgeAtMaturity (|r| > 0.7):\n")
print(high_cor_vars)

species_shift <- species_distr_stats %>%
  filter(Year >= 2080) %>%
  group_by(SpeciesPool, SpeciesID, Year) %>%
  summarise(
    AvgDiff = mean(diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Shift = ifelse(AvgDiff > 0, "Upward", "Downward"))

rf_shift <- niches %>%
  dplyr::select(-all_of(problem_cols)) %>%
  # - RecruitmentInvestmentRel, AgeAtMaturity
  #dplyr::select(-MaxRecruitsAtMaxMass, - MassAtMaturity) %>%  # Remove cols with similar contribution
  mutate(
    RangeTemp = MaxTemp - MinTemp,
    RangeHum = MaxHum - MinHum,
    RangeLight = MaxLight - MinLight
  ) %>%
  right_join(., species_shift, by = c("SpeciesPool", "SpeciesID")) %>%
  dplyr::select(-Shift) %>%
  drop_na() #%>% select(-SpeciesID, -Shift)

## set the seed to make your partition reproducible
set.seed(123)

train_idx <- sample(seq_len(nrow(rf_shift)), size = floor(0.75 * nrow(rf_shift)))

train <- rf_shift[train_idx,]
test <- rf_shift[-train_idx,]

rf_model <- randomForest(
  AvgDiff ~ ., data = train,
  importance = TRUE, ntree = 1000,
)

# --- 3. Model performance ----------------------------------------------

print(rf_model)
# Look at out-of-bag error rate and confusion matrix

# Get predicted values (continuous)
preds <- predict(rf_model, test)

# Compute residuals
residuals <- test$AvgDiff - preds

# Compute error metrics
MAE <- mean(abs(residuals))
RMSE <- sqrt(mean(residuals^2))
R2 <- 1 - (sum(residuals^2) / sum((test$AvgDiff - mean(test$AvgDiff))^2))

# Print results
cat("Mean Absolute Error (MAE):", round(MAE, 4), "\n")
cat("Root Mean Squared Error (RMSE):", round(RMSE, 4), "\n")
# Mean Absolute Error (MAE): 2.0839
# Root Mean Squared Error (RMSE): 2.9995

# --- 4. Variable importance --------------------------------------------

importance_df <- as.data.frame(importance(rf_model))
importance_df <- importance_df %>%
  rownames_to_column("Trait") %>%
  arrange(desc(`%IncMSE`))

# Print ranked trait importance
print(importance_df)

# --- 5. Visualize variable importance ----------------------------------

# Fix trait names
trait_labels <- c(
  "Year" = "Year",
  "MaxTemp" = "Max. temperature",
  "SpeciesID" = "Species identity",
  "OptimumTemp" = "Optimum remperature",
  "MinLight" = "Min. light",
  "GrowthRate" = "Growth rate",
  "RangeHum" = "Humidity range",
  "OptimumHum" = "Optimum humidity",
  "AgeAtMaturity" = "Age at maturity",
  "MinHum" = "Min. humidity",
  "DispersalKernel" = "Dispersal kernel",
  "RangeLight" = "Light range",
  "MaximumMass" = "Max. mass",
  "MaxHum" = "Max. humidity",
  "MassAtMaturity" = "Mass at maturity",
  "MaxLight" = "Max. light",
  "RecruitmentInvestmentRel" = "Rel. recruitment investment",
  "MaxRecruitsAtMaxMass" = "Max. recruits at max. mass",
  "SpeciesPool" = "Species pool",
  "MaxRecruitsAtMassAtMaturity" = "Max. recruits at maturity mass",
  "DispersalKernelAsymmetry" = "Dispersal kernel asymmetry",
  "OptimumLight" = "Optimum light",
  "RangeTemp" = "Temperature range",
  "MinTemp" = "Min. temperature"
)

# Apply to your plot
importance_df$Trait_display <- trait_labels[importance_df$Trait]
importance_df_surv$Trait_display <- trait_labels[importance_df_surv$Trait]
importance_df_surv <- importance_df_surv %>% drop_na()

# ----- Species shift analysis ------------------------------- #

# Plot importance
p1 <- ggplot(importance_df, aes(
    x = reorder(Trait_display, `%IncMSE`),
    y = `%IncMSE`
  )) +
  geom_col(fill = "#93A8AC") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_flip() +
  labs(
    x = "",
    y = "% increase in MSE (importance)"
  ) +
  theme_minimal() +     # <- No arguments here
  theme(
    text = element_text(size = 35),        # Axis labels, etc.
    axis.text.x = element_text(size = 30), # X-tick labels slightly smaller
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())

pdf(file.path(DirectoryPlots, "randomForestΤraitImportanceShiftTraits_v2.pdf"),
    width = 9, height = 10)
print(p1)
dev.off()

# Scatter plot Max. temperature vs. Avg. shift
pdf(file.path(DirectoryPlots, "scatter_shift_vs_maxtemp_v2.pdf"), width = 9, height = 8)
ggplot(rf_shift, aes(x = MaxTemp, y = AvgDiff)) +
  geom_point(color = "#93A8AC", size = 3, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_smooth(method = "lm", color = "#393838", se = TRUE) +
  labs(
    x = "Max. temperature (°C)",
    y = "Species shift with CC (m)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    axis.text.x = element_text(size = 25)
  )
dev.off()

# Compute correlation
cor.test(rf_shift$MaxTemp, rf_shift$AvgDiff, method = "pearson")
cor.test(rf_shift$MinLight, rf_shift$AvgDiff, method = "pearson")
cor.test(rf_shift$MaxHum, rf_shift$AvgDiff, method = "pearson")
cor.test(rf_shift$MaxLight, rf_shift$AvgDiff, method = "pearson")

# ----- Species survival analysis ------------------------------- #

p2 <- ggplot(importance_df_surv, aes(x = reorder(Trait_display, MeanDecreaseGini),
                          y = MeanDecreaseGini)) +
  geom_col(fill = "#93A8AC") +
  coord_flip() +
  labs(x = "",
       y = "Mean decrease\nin Gini (importance)") +
  theme_minimal() +
  theme(
    text = element_text(size = 35),        # Axis labels, etc.
    axis.text.x = element_text(size = 30), # X-tick labels slightly smaller
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())

pdf(file.path(DirectoryPlots, "randomForestΤraitImportanceSurvival_v1.pdf"), width = 10, height = 10)
print(p2)
dev.off()

# -- Max temp important
colors <- c('Survived CC' = '#C5D2D5', 'Died with CC' = '#E3E9EA')

p1 <- ggplot(rf_data, aes(x = survival, y = MaxTemp, fill = survival)) +
  geom_boxplot(notch = TRUE) +
  labs(x = "", y = "Max. temperature (°C)") +
  scale_x_discrete(labels = c(
    "Survived CC" = "Survived\nCC",
    "Died with CC" = "Died\nwith CC"
  )) +
  scale_fill_manual(
    name = "",
    values = c(
      "Survived CC" = colors[["Survived CC"]],
      "Died with CC" = colors[["Died with CC"]]
    )
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    legend.position = "none"
  )

p2 <- ggplot(rf_data, aes(x = survival, y = MaxHum, fill = survival)) +
  geom_boxplot(notch = TRUE) +
  labs(x = "", y = "Max. humidity (%)") +
  scale_x_discrete(labels = c(
    "Survived CC" = "Survived\nCC",
    "Died with CC" = "Died\nwith CC"
  )) +
  scale_fill_manual(
    name = "",
    values = c(
      "Survived CC" = colors[["Survived CC"]],
      "Died with CC" = colors[["Died with CC"]]
    )
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    legend.position = "none"
  )

# Save to PDF
pdf(file.path(DirectoryPlots, "boxplot_max_temp_hum_by_survival_CC_v4.pdf"),
    width = 9, height = 8)

(p1 | p2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 25)
  )

dev.off()

# Just because
p_light <- ggplot(rf_data, aes(x = survival, y = MinLight, fill = survival)) +
  geom_boxplot(notch = TRUE) +
  labs(x = "", y = "Min. light") +
  scale_x_discrete(labels = c(
    "Survived CC" = "Survived\nCC",
    "Died with CC" = "Died\nwith CC"
  )) +
  scale_fill_manual(
    name = "",
    values = c(
      "Survived CC" = colors[["Survived CC"]],
      "Died with CC" = colors[["Died with CC"]]
    )
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    legend.position = "none"
  )
pdf(file.path(DirectoryPlots, "boxplot_min_light_by_survival_CC_v4.pdf"), width = 5, height = 8)
print(p_light)
dev.off()