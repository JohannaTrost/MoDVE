library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)

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