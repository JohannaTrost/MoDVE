library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(patchwork) # for combining plots

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

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
  ) %>%
  drop_na(.)

summary_shift <- species_distr_stats %>%
  group_by(SpeciesPool, SpeciesID) %>%
  arrange(Year, .by_group = TRUE) %>%  # ensure data is sorted by Year
  reframe(
    diff = mean(tail(diff[!is.na(diff)], 20), na.rm = TRUE),
    CC = mean(tail(CC[!is.na(CC)], 20), na.rm = TRUE),
    last_year_alive = max(Year[!is.na(diff)]),
  ) %>%
  filter(last_year_alive >= 2050) %>%
  mutate(Shift = ifelse(diff > 0, "Upward", "Downward"))


# ---- Extract env. conditions

# - 1. Get time and location of individuals od interest
filtered_species_distr <- summary_shift %>%
  select(SpeciesPool, SpeciesID, last_year_alive, Shift, diff) %>%
  right_join(species_distr, ., by = c("SpeciesPool", "SpeciesID")) %>%
  filter(Year <= last_year_alive & Year >= last_year_alive - 20) %>%
  drop_na(.) # Drop all species that only survived in the sim. with community

# - 2. Loop over years and extract MC matrices and env. conditions
scenario_dir <- list("CC" = "climdata_era5_cmip6_1981-2100_ssp245",
                     "No CC" = "climdata_era5_cmip6_1981-2100_ssp245_no_cc")
species_env <- NA
for (scenario in c("CC", "No CC")) {
  for (forest in 0:2) {
    species_forest <- filtered_species_distr %>% filter(ForestID == forest, Scenario == scenario)

    for (ts in sort(unique(species_forest$TimeStep))) {
      species_locations <- species_forest %>% filter(TimeStep == ts)

      # Load mC matrix
      mh <- readRDS(file.path(base_dir, scenario_dir[[scenario]], "a1_2",
                              paste0("forest", forest), paste0("microhabitatMatrix", ts, ".rds")))
      mc <- mh[,,, c(3,4,5)]

      # Get the env conditions at the species position
      species_locations[, c("Light", "Hum", "Temp")] <- t(mapply(function(x, y, z) {
        mc[x, y, z, ]
      }, species_locations$X, species_locations$Y, species_locations$Z))

      # Get the minimum potential env. conditions on the vertical axes
      species_locations[, c("MinPotLight", "MinPotHum", "MinPotTemp")] <- t(mapply(function(x, y) {
        apply(mc[x, y, , ], 2, min)
      }, species_locations$X, species_locations$Y))

      # Get the max potential env. conditions on the vertical axes
      species_locations[, c("MaxPotLight", "MaxPotHum", "MaxPotTemp")] <- t(mapply(function(x, y) {
        apply(mc[x, y, , ], 2, max)
      }, species_locations$X, species_locations$Y))

      # Specifiy scenario
      species_locations$Scenario <- scenario

      if (all(is.na(species_env))) {
        species_env <- species_locations
      } else {
        species_env <- rbind(species_locations, species_env)
      }
    }
  }
}

# Now get the theoretical niche specifications
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

indiv_env_niches <- left_join(species_env, niches, by = c("SpeciesID", "SpeciesPool"))
#write_csv(indiv_env_niches, file.path(DirectoryPlots, "a5_individuals_cc_vs_no_cc_niches_envCond_shift.csv"))

# Species env. conditions in simulation with community
species_env_agg <- species_env %>%
  group_by(Scenario, SpeciesPool, SpeciesID, Shift) %>%
  summarize(AvgTempCond = mean(Temp, na.rm = TRUE),
            MaxTempCond = max(Temp, na.rm = TRUE),
            MinTempCond = min(Temp, na.rm = TRUE),
            AvgHumCond = mean(Hum, na.rm = TRUE),
            MaxHumCond = max(Hum, na.rm = TRUE),
            MinHumCond = min(Hum, na.rm = TRUE),
            AvgLightCond = mean(Light, na.rm = TRUE) * 900,
            MaxLightCond = max(Light, na.rm = TRUE) * 900,
            MinLightCond = min(Light, na.rm = TRUE) * 900,
            # Potential env. conditions per species
            MaxPotLight = mean(MaxPotLight, na.rm = TRUE) * 900,
            MaxPotHum = mean(MaxPotHum, na.rm = TRUE),
            MaxPotTemp = mean(MaxPotTemp, na.rm = TRUE),
            MinPotLight = mean(MinPotLight, na.rm = TRUE) * 900,
            MinPotHum = mean(MinPotHum, na.rm = TRUE),
            MinPotTemp = mean(MinPotTemp, na.rm = TRUE),
            AvgDiff = mean(diff, na.rm = TRUE),
            .groups = "drop")

# Merge realized and fundamental niches
species_env_niches <- left_join(species_env_agg, niches, by = c("SpeciesID", "SpeciesPool"))

# ---- Plot realised vs theoretical niches

# Assuming your tibble is named species_env_niches
df <- species_env_niches %>%
    filter(AvgDiff <= quantile(AvgDiff, 0.25) | AvgDiff >= quantile(AvgDiff, 0.75))

# --- Prepare data in long format for easier plotting ---
# We'll gather realised and theoretical values for Temp, Hum, Light

realised <- df %>%
  select(SpeciesPool, SpeciesID,
         AvgTempCond, MinTempCond, MaxTempCond,
         AvgHumCond, MinHumCond, MaxHumCond, Scenario,
         AvgLightCond, MinLightCond, MaxLightCond) %>%
  pivot_longer(
    cols = -c(SpeciesPool, SpeciesID, Scenario),
    names_to = c("Stat", "Variable"),
    names_pattern = "(Avg|Min|Max)(Temp|Hum|Light)Cond",
    values_to = "Value"
  ) %>%
  mutate(Type = ifelse(Scenario == "CC", "Realized niche (CC)", "Realized niche (no CC)"))

potential <- df %>%
  select(SpeciesPool, SpeciesID, MaxPotLight, MinPotLight, MaxPotHum, MinPotHum, MaxPotTemp,
         MinPotTemp, Scenario) %>%
  pivot_longer(
    cols = -c(SpeciesPool, SpeciesID, Scenario),
    names_to = c("Stat", "Variable"),
    names_pattern = "(Min|Max)Pot(Temp|Hum|Light)",
    values_to = "Value"
  ) %>%
  mutate(Type = ifelse(Scenario == "CC", "Env. range (CC)", "Env. range (no CC)"))

theoretical <- df %>%
  dplyr::select(SpeciesPool, SpeciesID, Scenario,
         OptimumTemp, MinTemp, MaxTemp,
         OptimumHum, MinHum, MaxHum,
         OptimumLight, MinLight, MaxLight) %>%
  pivot_longer(
    cols = -c(SpeciesPool, SpeciesID, Scenario),
    names_to = c("Stat", "Variable"),
    names_pattern = "(Optimum|Min|Max)(Temp|Hum|Light)",
    values_to = "Value"
  ) %>%
  mutate(Type = "Fundamental niche")

# Combine both
df_long <- bind_rows(realised, theoretical, potential) %>%
  mutate(
    Stat = recode(Stat,
                  "Avg" = "Mean", "Optimum" = "Mean"),
    Variable = recode(Variable,
                      "Temp" = "Temperature (°C)",
                      "Hum" = "Humidity (%)",
                      "Light" = "Light"),
    Label = paste0(SpeciesPool, "/", SpeciesID)
  )
df_wide <- df_long %>%
  pivot_wider(names_from = Stat, values_from = Value) %>%
  left_join(., df %>% select(SpeciesPool, SpeciesID, Shift, Scenario), by = c("SpeciesPool", "SpeciesID", "Scenario"))

# -- Plot realized vs fundamental niche
DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/position_shift")

pd <- position_dodge(width = 0.8)

# Control order of niche types shown
df_wide$Type <- factor(
  df_wide$Type,
  levels = c("Realized niche (CC)", "Env. range (CC)", "Realized niche (no CC)",
             "Env. range (no CC)", "Fundamental niche")
)

shift <- "Upward"
#shift <- "Downward"

# Filter to the current shift
df_shift <- df_wide %>% filter(Shift == shift)

# Determine order of species (Label) by max temperature for realized niche
order_labels <- df_shift %>%
  filter(Type == "Realized niche (CC)", Variable == "Temperature (°C)") %>%
  arrange(Max) %>%
  mutate(Label = factor(Label, levels = unique(Label))) %>%
  select(Label)

# Join the order back to your data
df_ordered <- df_shift %>%
  left_join(order_labels, by = c("Label")) %>%
  mutate(Label = factor(Label, levels = levels(order_labels$Label)))

#pdf(file.path(DirectoryPlots, paste0(shift, "_realised_vs_theoretical_vs_potential_niche_not_species.pdf")),
#    width = 18, height = 8)
ggplot(df_ordered,
       aes(x = Type, y = Mean, group = Label)) +
  geom_point(position = pd, size = 2) +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, alpha = 0.7, position = pd) +
  facet_wrap(~Variable, scales = "free_y", nrow = 3) +
  labs(
    title = "Temperature (°C) by SpeciesPool/SpeciesID",
    x = "SpeciesPool / SpeciesID",
    y = "Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
#dev.off()


# --- Compute range filling

range_filling <- species_env_niches %>%
  mutate(TempNicheFill = ((MaxTempCond - MinTempCond) / (MaxTemp - MinTemp)) * 100,
         HumNicheFill = ((MaxHumCond - MinHumCond) / (MaxHum - MinHum)) * 100,
         LightNicheFill = ((MaxLightCond - MinLightCond) / (MaxLight - MinLight)) * 100,
         TempRangeFill = ((MaxTempCond - MinTempCond) / (MaxPotTemp - MinPotTemp)) * 100,
         HumRangeFill = ((MaxHumCond - MinHumCond) / (MaxPotHum - MinPotHum)) * 100,
         LightRangeFill = ((MaxLightCond - MinLightCond) / (MaxPotLight - MinPotLight)) * 100)

# - Plot niche filling

# Reshape data to long format for ggplot
range_filling_long <- range_filling %>%
  select(Scenario, TempNicheFill, HumNicheFill, LightNicheFill) %>%
  pivot_longer(
    cols = c(TempNicheFill, HumNicheFill, LightNicheFill),
    names_to = "Variable",
    values_to = "NicheFill"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      TempNicheFill = "Temperature",
                      HumNicheFill = "Relative humidity",
                      LightNicheFill = "Light")
  )

# Define colors for scenarios
colors <- c('CC' = '#f7766e', 'No CC' = '#004aad')

write_csv(range_filling_long, file.path(base_dir, "a5_niche_filling_cc_vs_no_cc.csv"))

# Create the boxplot
output_file2 <- file.path(DirectoryPlots, "boxplot_niche_filling.pdf")
pdf(output_file2, width = 7, height = 6)
ggplot(range_filling_long, aes(x = Variable, y = NicheFill, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  labs(
    x = "",
    y = "Niche filling (%)",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()

# - Plot range filling

# Reshape data to long format for ggplot
range_filling_long <- range_filling %>%
  select(Scenario, TempRangeFill, HumRangeFill, LightRangeFill) %>%
  pivot_longer(
    cols = c(TempRangeFill, HumRangeFill, LightRangeFill),
    names_to = "Variable",
    values_to = "RangeFill"
  ) %>%
  mutate(
    Variable = recode(Variable,
                      TempRangeFill = "Temperature",
                      HumRangeFill = "Relative humidity",
                      LightRangeFill = "Light")
  )

# Define colors for scenarios
colors <- c('CC' = '#f7766e', 'No CC' = '#004aad')

# Create the boxplot
output_file2 <- file.path(DirectoryPlots, "boxplot_range_filling.pdf")
pdf(output_file2, width = 7, height = 6)
ggplot(range_filling_long, aes(x = Variable, y = RangeFill, fill = Scenario)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, notch = TRUE) +
  scale_fill_manual(values = colors) +
  labs(
    x = "",
    y = "Range filling (%)",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()

# ---- Relating shift and niche filling differences ---- #

# - Niche filling difference vs. shift

range_diff <- range_filling %>%
  select(SpeciesID, SpeciesPool, Scenario, AvgDiff, TempNicheFill, TempRangeFill,
         HumNicheFill, HumRangeFill, LightNicheFill, LightRangeFill, MaxTemp, MaxHum, MinLight) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = c(TempNicheFill, TempRangeFill, HumNicheFill, HumRangeFill, LightNicheFill,
                    LightRangeFill, MaxTemp, MaxHum, MinLight),
    names_sep = "_"
  ) %>%
  mutate(
    # Compute differences (CC - NoCC)
    Diff_TempNicheFill = TempNicheFill_CC - `TempNicheFill_No CC`,
    Diff_TempRangeFill = TempRangeFill_CC - `TempRangeFill_No CC`,
    Diff_HumNicheFill = HumNicheFill_CC - `HumNicheFill_No CC`,
    Diff_HumRangeFill = HumRangeFill_CC - `HumRangeFill_No CC`,
    Diff_LightNicheFill = LightNicheFill_CC - `LightNicheFill_No CC`,
    Diff_LightRangeFill = LightRangeFill_CC - `LightRangeFill_No CC`,
  )

# Compute correlartions
cor.test(range_diff$Diff_TempNicheFill, range_diff$AvgDiff, use = "complete.obs")
cor.test(range_diff$Diff_HumNicheFill, range_diff$AvgDiff, use = "complete.obs")
cor.test(range_diff$Diff_LightNicheFill, range_diff$AvgDiff, use = "complete.obs")

cor.test(range_diff$MaxTemp_CC, range_diff$AvgDiff, use = "complete.obs")

# Step 3: Scatter plot of shift vs niche filling difference
output_file2 <- file.path(DirectoryPlots, "scatter_shift_vs_niche_fill_diff_temp.pdf")

pdf(output_file2, width = 7, height = 6)
ggplot(range_diff, aes(x = Diff_TempNicheFill, y = AvgDiff)) +
  geom_point(size = 2, alpha = 0.7, color = "#004aad") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_smooth(method = "lm", se = TRUE, color = "#f7766e") +
  labs(
    x = expression(Delta~"Temperature niche filling (%)"~(CC - No~CC)),
    y = expression(Delta~"Position shift (m)"~(CC - No~CC))
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()
# Hum
output_file2 <- file.path(DirectoryPlots, "scatter_shift_vs_niche_fill_diff_hum.pdf")

pdf(output_file2, width = 7, height = 6)
ggplot(range_diff, aes(x = Diff_HumNicheFill, y = AvgDiff)) +
  geom_point(size = 2, alpha = 0.7, color = "#004aad") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_smooth(method = "lm", se = TRUE, color = "#f7766e") +
  labs(
    x = expression(Delta~"Relatve humidity niche filling (%)"~(CC - No~CC)),
    y = expression(Delta~"Position shift (m)"~(CC - No~CC))
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()
# Light
output_file2 <- file.path(DirectoryPlots, "scatter_shift_vs_niche_fill_diff_light.pdf")

pdf(output_file2, width = 7, height = 6)
ggplot(range_diff, aes(x = Diff_LightNicheFill, y = AvgDiff)) +
  geom_point(size = 2, alpha = 0.7, color = "#004aad") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_smooth(method = "lm", se = TRUE, color = "#f7766e") +
  labs(
    x = expression(Delta~"Light niche filling (%)"~(CC - No~CC)),
    y = expression(Delta~"Position shift (m)"~(CC - No~CC))
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )
dev.off()