library(lubridate)
library(readxl)
library(tidyverse)
library(lubridate)
library(readxl)
library(tidyverse)
library(ggplot2)


# Function to load and process a logger file
load_logger_data <- function(dir_name) {
  file_path <- file.path(emp_path, dir_name, "mc_data.xlsx")

  read_excel(file_path) %>%
    select("Date-Time (Brazil Standard Time)", "Temperature (°C)", "RH (%)") %>%
    rename(obs_time = "Date-Time (Brazil Standard Time)",
           tair = "Temperature (°C)",
           relhum = "RH (%)") %>%
    mutate(
      obs_time = force_tz(as.POSIXct(obs_time), tzone = "America/Sao_Paulo"),
      obs_time_utc = with_tz(obs_time, tzone = "UTC"),
      obs_hour_utc = ceiling_date(obs_time_utc, unit = "hour")
    ) %>%
    group_by(obs_hour_utc) %>%
    summarise(
      tair = mean(tair, na.rm = TRUE),
      relhum = mean(relhum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(obs_time = obs_hour_utc)
}


# Function to compare one logger against model and macroclimate at specific height
process_comparison_height <- function(dir_name, height_idx, mc_data_height) {
  emp_data <- load_logger_data(dir_name) %>%
    rename(tair_emp = tair, relhum_emp = relhum)

  # Join with model and macro
  joined <- emp_data %>%
    left_join(mc_data_height, by = "obs_time") %>%
    left_join(macro_data, by = "obs_time")

  tibble(
    logger = dir_name,
    height_idx = height_idx,

    # Model vs Empirical
    mae_tair_model = mean(abs(joined$tair_emp - joined$tair), na.rm = TRUE),
    cor_tair_model = cor(joined$tair_emp, joined$tair, use = "complete.obs"),
    mae_relhum_model = mean(abs(joined$relhum_emp - joined$relhum), na.rm = TRUE),
    cor_relhum_model = cor(joined$relhum_emp, joined$relhum, use = "complete.obs"),

    # Macroclimate vs Empirical (skip if same as macro reference)
    mae_tair_macro = if (dir_name != macro_dir) mean(abs(joined$tair_emp - joined$tair_macro), na.rm = TRUE) else NA_real_,
    cor_tair_macro = if (dir_name != macro_dir) cor(joined$tair_emp, joined$tair_macro, use = "complete.obs") else NA_real_,
    mae_relhum_macro = if (dir_name != macro_dir) mean(abs(joined$relhum_emp - joined$relhum_macro), na.rm = TRUE) else NA_real_,
    cor_relhum_macro = if (dir_name != macro_dir) cor(joined$relhum_emp, joined$relhum_macro, use = "complete.obs") else NA_real_
  )
}

# Function to find optimal height for all loggers
find_optimal_height <- function(emp_dirs, mc, obs_time, heights) {
  all_results <- tibble()

  # Test each height
  for (h in seq_along(heights)) {
    cat("Testing height", h, "of", length(heights), "\n")

    # Extract data for this height
    median_tair <- apply(mc$tair, c(3, 4), median, na.rm = TRUE)
    median_relhum <- apply(mc$relhum, c(3, 4), median, na.rm = TRUE)

    mc_data_height <- data.frame(
      obs_time = obs_time,
      tair = median_tair[, h],
      relhum = median_relhum[, h]
    )

    # Compare all loggers at this height
    height_results <- map_dfr(emp_dirs, ~process_comparison_height(.x, h, mc_data_height))
    all_results <- bind_rows(all_results, height_results)
  }

  return(all_results)
}

# Function to create vertical profile plot
create_vertical_profile <- function(mc, obs_time, emp_dirs, optimal_height, heights) {
  # Calculate average profiles across time for each height
  profile_data <- tibble()

  for (h in seq_along(heights)) {
    median_tair <- apply(mc$tair, c(3, 4), median, na.rm = TRUE)
    median_relhum <- apply(mc$relhum, c(3, 4), median, na.rm = TRUE)

    avg_tair <- mean(median_tair[, h], na.rm = TRUE)
    avg_relhum <- mean(median_relhum[, h], na.rm = TRUE)

    profile_data <- bind_rows(profile_data,
                              tibble(height = heights[h],
                                     tair_avg = avg_tair,
                                     relhum_avg = avg_relhum))
  }

  # Load empirical data for vertical lines
  emp_values <- tibble()
  for (dir in emp_dirs) {
    emp_data <- load_logger_data(dir) %>%
      summarise(
        logger = dir,
        tair_mean = mean(tair, na.rm = TRUE),
        relhum_mean = mean(relhum, na.rm = TRUE)
      )
    emp_values <- bind_rows(emp_values, emp_data)
  }

  # Create temperature profile plot
  p_temp <- ggplot(profile_data, aes(x = tair_avg, y = height)) +
    geom_line(color = "blue", size = 1.2) +
    geom_point(color = "blue", size = 2) +
    geom_vline(data = emp_values, aes(xintercept = tair_mean, color = logger),
               linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = heights[optimal_height], color = "red",
               linetype = "solid", size = 1, alpha = 0.8) +
    labs(x = "Temperature (°C)", y = "Height (m)",
         title = "Vertical Temperature Profile",
         subtitle = paste("Red line shows optimal height:", heights[optimal_height], "m")) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Logger Data"))

  # Create humidity profile plot
  p_humid <- ggplot(profile_data, aes(x = relhum_avg, y = height)) +
    geom_line(color = "green", size = 1.2) +
    geom_point(color = "green", size = 2) +
    geom_vline(data = emp_values, aes(xintercept = relhum_mean, color = logger),
               linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = heights[optimal_height], color = "red",
               linetype = "solid", size = 1, alpha = 0.8) +
    labs(x = "Relative Humidity (%)", y = "Height (m)",
         title = "Vertical Humidity Profile",
         subtitle = paste("Red line shows optimal height:", heights[optimal_height], "m")) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Logger Data"))

  return(list(temp_plot = p_temp, humid_plot = p_humid, profile_data = profile_data))
}

# --- Main Analysis ---

# Define paths
emp_path <- "/Users/johanna/Uni/masterarbeit/data/empirical/Datalogger 400m elevation REGUA understory Trilha Verde"

# Define logger directories
emp_dirs <- c("3600m, 446mASL", "3800m, 387mASL", "3650m, 438mASL", "3700m, 433mASL")
macro_dir <- "3850m, 387mASL waterfall reference"  # This is the macroclimate reference

# Load macroclimate data
macro_data <- load_logger_data(macro_dir) %>%
  rename(tair_macro = tair, relhum_macro = relhum)

# --- Load simulated data
start <- "2024-10-27"
stop <- "2024-10-31"
mc_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_output/regua_2024_test_v4"
mc <- readRDS(paste0(mc_dir, "/sim_with_era5_", start, "_", stop, ".rds"))

# Create time sequence
obs_time <- seq(ymd_hms(paste(start, "00:00:00")), ymd_hms(paste(stop, "23:59:59")), by = "hour")

# Define height levels (adjust based on your model structure)
# Assuming heights are in meters - modify these values based on your actual height levels
heights <- seq(0.5, dim(mc$tair)[4], by = 1)  # Example: from 0.1m to 10m in 0.5m steps
# Or if you know the exact heights: heights <- c(0.1, 0.5, 1, 2, 3, 5, 10)

# Find optimal height across all loggers
cat("Finding optimal height across all loggers...\n")
all_height_results <- find_optimal_height(emp_dirs, mc, obs_time, heights)

# Calculate overall performance metrics for each height
height_summary <- all_height_results %>%
  group_by(height_idx) %>%
  summarise(
    height = heights[height_idx[1]],
    mean_mae_tair = mean(mae_tair_model, na.rm = TRUE),
    mean_mae_relhum = mean(mae_relhum_model, na.rm = TRUE),
    mean_cor_tair = mean(cor_tair_model, na.rm = TRUE),
    mean_cor_relhum = mean(cor_relhum_model, na.rm = TRUE),
    combined_error = mean_mae_tair + mean_mae_relhum,  # Combined error metric
    .groups = 'drop'
  ) %>%
  arrange(combined_error)

# Select optimal height (minimum combined error)
optimal_height_idx <- height_summary$height_idx[1]
optimal_height <- heights[optimal_height_idx]

cat("Optimal height found:", optimal_height, "m\n")
cat("Combined MAE at optimal height:", height_summary$combined_error[1], "\n")

# Show height summary
print("Height Performance Summary (top 10):")
print(head(height_summary, 10))

# Get results for optimal height
optimal_results <- all_height_results %>%
  filter(height_idx == optimal_height_idx) %>%
  arrange(mae_relhum_model)

print("\nResults at optimal height:")
print(optimal_results)

# Create vertical profile plots
cat("Creating vertical profile plots...\n")
profile_plots <- create_vertical_profile(mc, obs_time, emp_dirs, optimal_height_idx, heights)

# Display plots
print(profile_plots$temp_plot)
print(profile_plots$humid_plot)

# Save plots if desired
ggsave("../../figs/mc_output/temperature_profile_mc_regua.png",
       profile_plots$temp_plot, width = 10, height = 8, dpi = 300)
ggsave("../../figs/mc_output/humidity_profile_mc_regua.png",
       profile_plots$humid_plot, width = 10, height = 8, dpi = 300)

# Create combined plot
library(gridExtra)
combined_plot <- grid.arrange(profile_plots$temp_plot, profile_plots$humid_plot, ncol = 2)

# Return key results
list(
  optimal_height = optimal_height,
  optimal_height_idx = optimal_height_idx,
  height_summary = height_summary,
  optimal_results = optimal_results,
  profile_data = profile_plots$profile_data,
  plots = profile_plots
)

# --- Compare empirical data against model and macroclimate

# Apply to all loggers
results <- map_dfr(emp_dirs, process_comparison)

# View sorted by model RH MAE
results_sorted <- results %>% arrange(mae_relhum_model)
print(results_sorted)

# --- Show results of the comparison (MAE and correlation)

# Summarize across loggers (excluding macro itself)
summary_stats <- results %>%
  filter(logger != macro_dir) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), na.rm = TRUE)) %>%
  t() %>%
  round(2)

print(summary_stats)

# --- Some plots

# Load and process all logger data
emp_data_list <- map(emp_dirs, load_logger_data)

# Combine all into a single tibble by full_join by time
repl_tair <- paste0("tair_", emp_dirs[length(emp_dirs)])
repl_relhum <- paste0("relhum_", emp_dirs[length(emp_dirs)])

emp_data_combined <- emp_data_list %>%
  reduce(full_join, by = "obs_time") %>%
  rename_at(vars(starts_with("tair")),
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) %>%
  rename_at(vars(starts_with("relhum")),
            str_replace, "(\\.\\w)+$", paste0("_", unlist(emp_dirs))) %>%
  rename("tair" := !!repl_tair, "relhum" := !!repl_relhum)

names(emp_data_combined)

# Compute mean and standard deviation across loggers for each time point
emp_summary <- emp_data_combined %>%
  rowwise() %>%
  mutate(
    tair_mean = mean(c_across(starts_with("tair_")), na.rm = TRUE),
    tair_sd = sd(c_across(starts_with("tair_")), na.rm = TRUE),
    relhum_mean = mean(c_across(starts_with("relhum_")), na.rm = TRUE),
    relhum_sd = sd(c_across(starts_with("relhum_")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(obs_time, tair_mean, tair_sd, relhum_mean, relhum_sd) %>%
  inner_join(., macro_data, by = "obs_time") %>%
  inner_join(., mc_sim, by = "obs_time")

# First plot: Air temperature (empirical vs. simulated)
plot_airt <- ggplot(emp_summary) +
  geom_ribbon(aes(x = obs_time,
                  ymin = tair_mean - 1.96 * tair_sd,
                  ymax = tair_mean + 1.96 * tair_sd),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = obs_time, y = tair_mean, color = "Mean Emp. MC"),
            size = 1.2) +
  geom_line(aes(x = obs_time, y = tair_macro, color = "Macroclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = tair, color = "Simulated"),
            size = 1) +
  scale_color_manual(name = "",
                     values = c("Mean Emp. MC" = "#2E86AB",     # Blue
                                "Macroclimate" = "#A23B72",     # Purple/magenta
                                "Simulated" = "#F18F01"),       # Orange
                     breaks = c("Mean Emp. MC", "Macroclimate", "Simulated")) +
  labs(title = "Air Temperature: Empirical vs. Simulated",
       x = "Observation Time",
       y = "Air Temperature (°C)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Second plot: Relative Humidity (empirical vs. simulated)
plot_relhum <- ggplot(emp_summary) +
  geom_ribbon(aes(x = obs_time,
                  ymin = relhum_mean - 1.96 * relhum_sd,
                  ymax = relhum_mean + 1.96 * relhum_sd),
              fill = "#2E86AB", alpha = 0.3) +
  geom_line(aes(x = obs_time, y = relhum_mean, color = "Mean Emp. MC"),
            size = 1.2) +
  geom_line(aes(x = obs_time, y = relhum_macro, color = "Macroclimate"),
            size = 1) +
  geom_line(aes(x = obs_time, y = relhum, color = "Simulated"),
            size = 1) +
  scale_color_manual(name = "",
                     values = c("Mean Emp. MC" = "#2E86AB",     # Blue
                                "Macroclimate" = "#A23B72",     # Purple/magenta
                                "Simulated" = "#F18F01"),       # Orange
                     breaks = c("Mean Emp. MC", "Macroclimate", "Simulated")) +
  labs(title = "Relative Humidity: Empirical vs. Simulated",
       x = "Observation Time",
       y = "Relative Humidity (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Print plots to a pdf file
pdf("../../figs/mc_output/v2_airt_emp_vs_sim_mc_regua.pdf")
print(plot_airt)
dev.off()

pdf("../../figs/mc_output/v2_relhum_emp_vs_sim_mc_regua.pdf")
print(plot_relhum)
dev.off()

pdf("../../figs/mc_output/v2_relhum_emp_vs_sim_mc_regua_ccf.pdf")
ccf(emp_sim_data$tair, emp_sim_data$tair_emp)
dev.off()
