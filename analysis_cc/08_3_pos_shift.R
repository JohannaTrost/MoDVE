library(readr)
library(dplyr)
library(tidyr)
library(ade4)
library(factoextra)
library(vegan)

base_dir <- file.path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")

species_distr <- read_csv(file.path(base_dir, "a5_species_distribution_cc_vs_no_cc.csv"))

species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

# Plot 2030 pos
DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/position_shift")

pdf(file.path(DirectoryPlots, "boxplot_position_2080-2100.pdf"), width = 4, height = 4)
species_distr_stats %>%
  filter(Year >= 2080) %>%
ggplot(., aes(x = Scenario, y = Position, fill = Scenario)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Scenario",
    y = "Species position (m)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
dev.off()

year <- 2080
species_shift <- species_distr_stats %>%
  filter(Year >= year) %>%
  group_by(SpeciesPool, SpeciesID) %>%
  summarise(
    AvgDiff = mean(diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Shift = ifelse(AvgDiff > 0, "Upward", "Downward")) %>%
  drop_na(.)


# Same plot without communitiy just upward shifting species

species_distr <- read_csv(file.path(base_dir, "a5_upward_shifted_species_distribution_cc_vs_no_cc.csv"))

species_distr_stats <- species_distr %>%
  group_by(Scenario, SpeciesPool, SpeciesID, ForestID, TimeStep, Year) %>%
  summarize(Position = mean(Height, na.rm = TRUE),
            Mass = mean(Mass, na.rm = TRUE),
            IQR = IQR(Height, na.rm = TRUE),
            Range = max(Height, na.rm = TRUE) - min(Height, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Scenario = factor(Scenario, levels = c("No CC", "CC"))) %>%
  arrange(Scenario)

# Plot 2030 pos
DirectoryPlots <- file.path("../../figs/a5_plots_test/cc_vs_no_cc/position_shift")

pdf(file.path(DirectoryPlots, "boxplot_position_no_comm_2080-2100.pdf"), width = 4, height = 4)
species_distr_stats %>%
  filter(Year >= 2080) %>%
ggplot(., aes(x = Scenario, y = Position, fill = Scenario)) +
  geom_boxplot(notch = TRUE) +
  labs(
    x = "Scenario",
    y = "Species position (m)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
dev.off()