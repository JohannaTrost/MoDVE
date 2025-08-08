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
