library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)
library(vegan)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

# Get species niches
niches <- NULL
for (sp in 1:10) {
  sp_niches <- read_csv(file.path(base_dir, "a2_1", paste0("SpeciesPool", sp, ".csv")),
                        show_col_types = FALSE) %>%
    dplyr:: select(-LightResponseA, -LightResponseB, -LightResponseC, -MinWind, -MaxWind, -OptimumWind,
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
  dplyr:: select(-SpeciesID, -SpeciesPool) %>%
  mutate(across(everything(), ~ scale(.)[,1])) %>%
  bind_cols(niches %>% dplyr:: select(SpeciesID, SpeciesPool), .)

# compute sds (NA if column is all NA)
sds <- sapply(niches_scaled, sd, na.rm = TRUE)
# columns with zero sd or NA sd
problem_cols <- names(sds)[is.na(sds) | sds == 0]

# PCA
pca_data <- niches_scaled %>%
  dplyr::dplyr:: select(-all_of(problem_cols)) %>%
  dplyr::dplyr:: select(-RecruitmentInvestmentRel, -MaxRecruitsAtMaxMass, -AgeAtMaturity, - MassAtMaturity) %>%  # Remove cols with similar contribution
  mutate(
    RangeTemp = MaxTemp - MinTemp,
    RangeHum = MaxHum - MinHum,
    RangeLight = MaxLight - MinLight
  )

res.pca <- pca_data %>%
  dplyr::dplyr:: select(-SpeciesID, -SpeciesPool) %>%
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
  dplyr::dplyr:: select(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year, Position) %>%
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
      !has_CC & has_NoCC ~ "Survived CC",
      has_CC & has_NoCC  ~ "Died with CC",
      TRUE ~ "none"  # optional: if both are NA
    )
  ) %>%
  dplyr:: select(SpeciesPool, SpeciesID, survival)

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
rf_data$survival <- as.factor(rf_data$survival)
trait_cols <- names(plt_data)[!(names(plt_data) %in%
  c("SpeciesID", "survival", "RecruitmentInc", "MaximmumMass"))]
rf_data$Survival <- ifelse(rf_data$survival %in% c("only CC", "both"), "Survived CC", "Died with CC")

# --- 2. RF model to predict survival --------------------------------
rf_model_balanced <- randomForest(
  survival ~ ., data = rf_data[, c("survival", trait_cols)],
  importance = TRUE, ntree = 1000,
  sampsize = rep(min(table(rf_data$survival)), 3)
)

# --- 3. Model performance ----------------------------------------------

print(rf_model_balanced)
# Look at out-of-bag error rate and confusion matrix

# --- 4. Variable importance --------------------------------------------

importance_df <- as.data.frame(importance(rf_model_balanced))
importance_df <- importance_df %>%
  rownames_to_column("Trait") %>%
  arrange(desc(MeanDecreaseGini))

# Print ranked trait importance
print(importance_df)

# --- 5. Visualize variable importance ----------------------------------

pdf(file.path(DirectoryPlots, "randomForestΤraitImportanceSurvival.pdf"), width = 7, height = 5)
ggplot(importance_df, aes(x = reorder(Trait, MeanDecreaseGini),
                          y = MeanDecreaseGini)) +
  geom_col(fill = "#56B4E9") +
  coord_flip() +
  labs(x = "Trait",
       y = "Mean Decrease in Gini (importance)") +
  theme_minimal(base_size = 13)
dev.off()

# --- 6. Optional: partial dependence (to interpret top traits) ----------

# For the most important traits:
library(pdp)

top_traits <- importance_df$Trait[1:5]
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

# -- Max temp important
for (trait in top_traits) {
  pdf(file.path(DirectoryPlots, paste0("boxplot_", trait, "_by_survival_CC_v2.pdf")),
      width = 6, height = 4)
  p <- ggplot(rf_data, aes(x = Survival, y = .data[[trait]], fill = Survival)) +
    geom_boxplot(notch = TRUE) +
    labs(
      x = "Survival",
      y = trait
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
  print(p)
  dev.off()
}

# ----------- Mortality

species_mortality_cat <- species_distr_stats %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    # Alive status between 2000–2020
    alive_NoCC_early = any(!is.na(`No CC`) & Year >= 2000 & Year <= 2020),
    alive_CC_early   = any(!is.na(CC) & Year >= 2000 & Year <= 2020),

    # Alive status from 2080 onward
    alive_NoCC_late = any(!is.na(`No CC`) & Year >= 2080),
    alive_CC_late   = any(!is.na(CC) & Year >= 2080),
    .groups = "drop"
  ) %>%
  mutate(
    # Mortality logic: alive before, but not anymore after 2080
    has_NoCC_mortality = alive_NoCC_early & !alive_NoCC_late,
    has_CC_mortality   = alive_CC_early & !alive_CC_late,

    mortality = case_when(
      has_CC_mortality & !has_NoCC_mortality ~ "only CC",
      !has_CC_mortality & has_NoCC_mortality ~ "only No CC",
      has_CC_mortality & has_NoCC_mortality  ~ "both",
      TRUE ~ "none"  # survived in both or never alive
    )
  ) %>%
  dplyr:: select(SpeciesPool, SpeciesID, mortality)

# Merge with niches
pca_data_mort <- pca_data %>%
  left_join(., species_mortality_cat, by = c("SpeciesPool", "SpeciesID")) %>%
  mutate(mortality = ifelse(is.na(mortality), "none", mortality))

# - Biplot with mortality categories

# Add a grouping column with NAs for excluded individuals
pca_data_mort$mortality_filtered <- pca_data_mort$mortality
pca_data_mort$mortality_filtered <- ifelse(
  pca_data_mort$mortality == "none",
  "Excluded",
  as.character(pca_data_mort$mortality)
)

pdf(file.path(DirectoryPlots, "trait_pca_biplot_mortality.pdf"))
fviz_pca_biplot(res.pca,
  label = "var",
  labelsize = 3,
  habillage = factor(pca_data_mort$mortality_filtered),
  palette = c("#00AFBB", "#FFFFFF", "#FC4E07", "#EFD2CB", "#241623"),
  addEllipses = TRUE,
  ellipse.type = "confidence",
  legend.title = "Mortality Groups",
  repel = TRUE
)
dev.off()

# -- Plot denisty only

# Extract PCA coordinates for individuals
pca_coords <- as.data.frame(res.pca$li[, 1:2])  # First two PCs
colnames(pca_coords) <- c("PC1", "PC2")

# Add survival information
pca_coords$Mortality <- pca_data_mort$mortality
pca_coords <- pca_coords %>% filter(Mortality != "none")

# Remove missing values (if any)
pca_coords <- na.omit(pca_coords)

# Compute group centroids
centroids <- pca_coords %>%
  group_by(Mortality) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE)
  )

# Compute explained variance (percent)
eig_var <- res.pca$eig / sum(res.pca$eig) * 100

pdf(file.path(DirectoryPlots, "trait_pca_mortality_group_centroids_density.pdf"), width = 7, height = 6)

ggplot(pca_coords, aes(x = PC1, y = PC2, color = Mortality, fill = Mortality)) +
  # Convex hulls
  geom_polygon(
    data = do.call(rbind, lapply(split(pca_coords, pca_coords$Mortality), function(df) {
      df[chull(df$PC1, df$PC2), ]
    })),
    aes(group = Mortality),
    alpha = 0.15,
    color = NA
  ) +
  # Density contours
  geom_density_2d(alpha = 0.3) +
  # Centroids
  geom_point(data = centroids, aes(x = PC1, y = PC2), shape = 21, size = 5, fill = "white", color = "black") +
  geom_text(data = centroids, aes(label = Mortality), vjust = -1, fontface = "bold") +
  # Axis labels with variance explained
  labs(
    x = paste0("PC1 (", round(eig_var[1], 1), "%)"),
    y = paste0("PC2 (", round(eig_var[2], 1), "%)"),
    color = "Mortality Groups",
    fill = "Mortality Groups"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

dev.off()

# Conclusion:
# Mortality patterns among species are largely overlapping across climate scenarios, indicating that
# major trait axes (growth rate, body size, humidity, and temperature tolerance) do not strongly
# distinguish which species die only under climate change versus those dying without it. However,
# there’s a hint that climate-change-specific mortality may be more associated with species having
# higher humidity dependence or narrower temperature tolerances.

# Join mortality PCS with categroy
pca_data_mort <- pca_data_mort %>%
  mutate(
    PC1 = res.pca$li$Axis1,
    PC2 = res.pca$li$Axis2,
    PC3 = res.pca$li$Axis3,
    PC4 = res.pca$li$Axis4,
    PC5 = res.pca$li$Axis5
  )


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

pdf(file.path(DirectoryPlots, paste0("trait_pca_biplot_shift_", year, "-2100.pdf")))
fviz_pca_biplot(res.pca,
                label = "var",
                labelsize = 3,
             dplyr:: select.ind = list(ind = na_idx),
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

