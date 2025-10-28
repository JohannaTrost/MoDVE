library(readr)
library(dplyr)
library(tidyr)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

species_distr <- read_csv(file.path(base_dir, "a5_upward_shifted_species_distribution_cc_vs_no_cc.csv"))

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
    initial_diff = mean(diff[Year == 2000], na.rm = TRUE),
    diff = mean(tail(diff[!is.na(diff)], 20), na.rm = TRUE),
    CC = mean(tail(CC[!is.na(CC)], 20), na.rm = TRUE),
    last_year_alive = max(Year[!is.na(diff)]),
  ) %>%
  filter(last_year_alive >= 2050) %>%
  mutate(Shift = ifelse(diff > 0, "Upward", "Downward"))

# --- Compare to shift with community

upward_shift_spec <- read_csv(file.path(base_dir, "a5_species_shift_upwards_cc_vs_no_cc.csv"))

shift_compare <- full_join(summary_shift, upward_shift_spec, by = c("SpeciesPool", "SpeciesID"))

DirectoryPlots <- file.path("../../../figs/a5_plots_test/cc_vs_no_cc/position_shift")
dir.create(DirectoryPlots)

# - Now comparing CC and no CC

lims <- range(c(shift_compare$diff, shift_compare$diff_comm), na.rm = TRUE)

pdf(file.path(DirectoryPlots, "pos_shift_with_cc_comm_vs_nocomm_upward_species.pdf"), width = 10, height = 8)
ggplot(shift_compare, aes(x = diff, color = Shift)) +
  geom_point(aes(y = diff_comm)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(
    x = "Species position shift without community",
    y = "Species position shift with community",
    color = "Shift"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Compute correlation between GrowthRate and shift for the final community species
cor_test <- cor.test(shift_compare$diff,
                     shift_compare$diff_comm,
                     method = "pearson")
cor_test # -> no correlation

# - Color by Optimal temp
pdf(file.path(DirectoryPlots, "pos_shift_with_cc_comm_vs_nocomm_upward_species_by_opthum.pdf"), width = 10, height = 8)
ggplot(shift_compare, aes(x = diff, color = OptimumHum)) +
  geom_point(aes(y = diff_comm)) +
  labs(
    x = "Species position shift without community",
    y = "Species position shift with community",
    color = "Optimum humidity"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()


# ---- Extract env. conditions

# - 1. Get time and location of individuals od interest
filtered_species_distr <- species_distr %>%
  right_join(.,
    shift_compare %>%
      select(SpeciesPool, SpeciesID, last_year_alive, Shift),
    by = c("SpeciesPool", "SpeciesID")
  ) %>%
  filter(Year <= last_year_alive & Year >= last_year_alive - 20, Shift == "Upward") %>%
  drop_na(.) # Drop all species that only survived in the sim. with community

# - 2. Loop over years and extract MC matrices and env. conditions
species_env <- NA
for (forest in 0:2) {
  species_forest <- filtered_species_distr %>% filter(ForestID == forest, Scenario == "CC")

  for (ts in sort(unique(species_forest$TimeStep))) {
    species_locations <- species_forest %>% filter(TimeStep == ts)

    # Load mC matrix
    mh <- readRDS(file.path(base_dir, "climdata_era5_cmip6_1981-2100_ssp245", "a1_2",
                            paste0("forest", forest), paste0("MicrohabitatMatrix", ts, ".rds")))
    mc <- mh[,,, c(3,4,5)]
    species_locations[, c("Light", "Hum", "Temp")] <- t(mapply(function(x, y, z) {
      mc[x, y, z, ]
    }, species_locations$X, species_locations$Y, species_locations$Z))

    if (all(is.na(species_env))) {
      species_env <- species_locations
    } else {
      species_env <- rbind(species_locations, species_env)
    }
  }
}

# Species env. conditions in simulation without community
species_env_agg <- species_env %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarize(AvgTempCond = mean(Temp, na.rm = TRUE),
            MaxTempCond = max(Temp, na.rm = TRUE),
            MinTempCond = min(Temp, na.rm = TRUE),
            AvgHumCond = mean(Hum, na.rm = TRUE),
            MaxHumCond = max(Hum, na.rm = TRUE),
            MinHumCond = min(Hum, na.rm = TRUE),
            AvgLightCond = mean(Light, na.rm = TRUE) * 900,
            MaxLightCond = max(Light, na.rm = TRUE) * 900,
            MinLightCond = min(Light, na.rm = TRUE) * 900,
            .groups = "drop")

# Now get the theoretical niche specifications
species_env_niches <- upward_shift_spec %>%
  select(MinTemp, MaxTemp, OptimumTemp,
         MinHum, MaxHum, OptimumHum,
         MinLight, MaxLight, OptimumLight,
         SpeciesPool, SpeciesID
  ) %>%
  left_join(species_env_agg, ., by = c("SpeciesID", "SpeciesPool"))

# ---- Plot realised vs theoretical niches

library(tidyverse)
library(patchwork) # for combining plots

# Assuming your tibble is named species_env_niches
df <- species_env_niches #%>%
  # mutate(MinLight = MinLight * 900,
  #        MaxLight = MaxLight * 900,
  #        OptimumLight = OptimumLight * 900
  # )

# --- Prepare data in long format for easier plotting ---
# We'll gather realised and theoretical values for Temp, Hum, Light

realised <- df %>%
  select(SpeciesPool, SpeciesID,
         AvgTempCond, MinTempCond, MaxTempCond,
         AvgHumCond, MinHumCond, MaxHumCond,
         AvgLightCond, MinLightCond, MaxLightCond) %>%
  pivot_longer(
    cols = -c(SpeciesPool, SpeciesID),
    names_to = c("Stat", "Variable"),
    names_pattern = "(Avg|Min|Max)(Temp|Hum|Light)Cond",
    values_to = "Value"
  ) %>%
  mutate(Type = "Realised")

theoretical <- df %>%
  select(SpeciesPool, SpeciesID,
         OptimumTemp, MinTemp, MaxTemp,
         OptimumHum, MinHum, MaxHum,
         OptimumLight, MinLight, MaxLight) %>%
  pivot_longer(
    cols = -c(SpeciesPool, SpeciesID),
    names_to = c("Stat", "Variable"),
    names_pattern = "(Optimum|Min|Max)(Temp|Hum|Light)",
    values_to = "Value"
  ) %>%
  mutate(Type = "Theoretical")

# Combine both
df_long <- bind_rows(realised, theoretical) %>%
  mutate(
    Stat = recode(Stat,
                  "Avg" = "Mean", "Optimum" = "Mean"),
    Variable = recode(Variable,
                      "Temp" = "Temperature (°C)",
                      "Hum" = "Humidity (%)",
                      "Light" = "Light"),
    Label = paste0(SpeciesPool, "/", SpeciesID)
  )
df_wide <- df_long %>% pivot_wider(names_from = Stat, values_from = Value)

# Display
pd <- position_dodge(width = 0.6)
pdf(file.path(DirectoryPlots, "no_comm_upward_realised_vs_theoretical_niche.pdf"))
ggplot(df_wide, aes(x = Label, y = Mean, color = Type, group = Type)) +
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
dev.off()
