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
    select(-LightResponseA, -LightResponseB, -LightResponseC, -MinWind, -MaxWind, -OptimumWind,
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
  select(-SpeciesID, -SpeciesPool) %>%
  mutate(across(everything(), ~ scale(.)[,1])) %>%
  bind_cols(niches %>% select(SpeciesID, SpeciesPool), .)

# compute sds (NA if column is all NA)
sds <- sapply(niches_scaled, sd, na.rm = TRUE)
# columns with zero sd or NA sd
problem_cols <- names(sds)[is.na(sds) | sds == 0]

# PCA
pca_data <- niches_scaled %>%
  select(-all_of(problem_cols)) %>%
  select(-RecruitmentInvestmentRel, -MaxRecruitsAtMaxMass, -AgeAtMaturity, - MassAtMaturity) %>%  # Remove cols with similar contribution
  mutate(
    RangeTemp = MaxTemp - MinTemp,
    RangeHum = MaxHum - MinHum,
    RangeLight = MaxLight - MinLight
  ) %>%
  select(-MaxTemp, -MinTemp, -MaxHum, -MinHum, -MaxLight, -MinLight)

res.pca <- pca_data %>%
  select(-SpeciesID, -SpeciesPool) %>%
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

pdf(file.path(DirectoryPlots, "trait_pca_biplot_vars.pdf"))
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
      has_CC & !has_NoCC ~ "only CC",
      !has_CC & has_NoCC ~ "only No CC",
      has_CC & has_NoCC  ~ "both",
      TRUE ~ "none"  # optional: if both are NA
    )
  ) %>%
  select(SpeciesPool, SpeciesID, survival)

# Merge with niches
pca_data_surv <- pca_data %>%
  left_join(., species_survival_cat, by = c("SpeciesPool", "SpeciesID")) %>%
  mutate(survival = ifelse(is.na(survival), "none", survival))

# Biplot with survival categories

# Add a grouping column with NAs for excluded individuals
survived_idx <- which(pca_data_surv$survival != "none")
pca_data_surv$survival_filtered <- pca_data_surv$survival
pca_data_surv$survival_filtered[-survived_idx] <- NA

pdf(file.path(DirectoryPlots, "trait_pca_biplot_survival.pdf"))
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

# Check if groups differ significantly in PCA space
pca_data_surv <- pca_data_surv %>%
  mutate(
    PC1 = res.pca$li$Axis1,
    PC2 = res.pca$li$Axis2,
    PC3 = res.pca$li$Axis3,
    PC4 = res.pca$li$Axis4,
    PC5 = res.pca$li$Axis5
  )

manova_test <- manova(cbind(PC1, PC2, PC3, PC4, PC5) ~ survival, data = pca_data_surv)
summary(manova_test)          # Pillai's trace test
summary.aov(manova_test)      # Separate univariate ANOVAs


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
  select(SpeciesPool, SpeciesID, mortality)

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

year <- 2030
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
