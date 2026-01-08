library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

mc_matrix <- readRDS("/Users/johanna/Uni/masterarbeit/data/mc_output/v3_mc_matrix.rds")

# Define variable names and units for better labeling
var_names <- c(
  "Mean Annual Temperature (°C)",
  "Mean Diurnal Temperature Range (°C)",
  "Annual Temperature Range (°C)",
  "Max Temperature of Warmest Month (°C)",
  "Min Temperature of Coldest Month (°C)",
  "Isothermality (%)",
  "Mean Annual Relative Humidity (%)",
  "Annual Humidity Range (%)",
  "Max Humidity (%)",
  "Min Humidity (%)",
  "Mean Annual Wind Speed (m/s)",
  "Annual Wind Speed Range (m/s)",
  "Max Wind Speed (m/s)",
  "Min Wind Speed (m/s)"
)

# Convert height index to actual height in meters
# Height 1 = 0.5m, then increases by 1m: 0.5, 1.5, 2.5, ..., 58.5
heights <- seq(0.5, 58.5, by = 1)

# Create a function to plot individual vertical profiles
plot_vertical_profile <- function(var_index, x_coord, y_coord) {
  # Extract data for this variable
  values <- mc_profile[x_coord, y_coord, 1:length(heights), var_index]

  # Create data frame
  profile_data <- data.frame(
    height = heights,
    value = values
  )

  # Create the plot
  ggplot(profile_data, aes(x = value, y = height)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 1.5, alpha = 0.7) +
    labs(
      title = var_names[var_index],
      x = "Value",
      y = "Height (m)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    ) +
    scale_y_continuous(breaks = seq(0, 60, by = 10))
}

# Create all 14 plots (assuming x and y coordinates - you'll need to specify these)
# For demonstration, I'll use x = 25, y = 25 (middle of your 50x50 grid)
x_coord <- 1  # Change this to your desired x coordinate
y_coord <- 1  # Change this to your desired y coordinate
n_temp_metrics <- 14
max_hgt <- 59

mc_profile <- array(rep(NA, 1 * 1 * 80 * n_temp_metrics),
                     dim = c(1, 1, 80, n_temp_metrics))
mc_profile[1, 1, , ] <- apply(mc_matrix, c(3, 4), mean, na.rm = TRUE)

# Generate all plots
plot_list <- list()
for (i in 1:14) {
  plot_list[[i]] <- plot_vertical_profile(i, x_coord, y_coord)
}

# Option 1: Display all plots in a grid (might be crowded)
grid_plot <- grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)

# Option 2: Create separate plots for different variable groups
# Temperature variables (1-6)
temp_plots <- grid.arrange(grobs = plot_list[1:6], ncol = 3, nrow = 2,
                          top = "Temperature Variables - Vertical Profiles")

# Humidity variables (7-10)
humidity_plots <- grid.arrange(grobs = plot_list[7:10], ncol = 2, nrow = 2,
                              top = "Humidity Variables - Vertical Profiles")

# Wind speed variables (11-14)
wind_plots <- grid.arrange(grobs = plot_list[11:14], ncol = 2, nrow = 2,
                          top = "Wind Speed Variables - Vertical Profiles")

# Option 3: Create a combined data frame for faceted plotting
create_combined_data <- function(x_coord, y_coord) {
  combined_data <- data.frame()

  for (i in 1:14) {
    values <- mc_profile[x_coord, y_coord, 1:max_hgt, i]
    temp_data <- data.frame(
      height = heights,
      value = values,
      variable = var_names[i],
      var_group = case_when(
        i <= 6 ~ "Temperature",
        i <= 10 ~ "Humidity",
        TRUE ~ "Wind Speed"
      )
    )
    combined_data <- rbind(combined_data, temp_data)
  }
  return(combined_data)
}

# Create faceted plot
combined_data <- create_combined_data(x_coord, y_coord)

faceted_plot <- ggplot(combined_data, aes(x = value, y = height)) +
  geom_line(color = "blue", size = 0.8) +
  geom_point(color = "red", size = 1, alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  labs(
    title = paste("Vertical Profiles of microclimate Variables (Grid Point:", x_coord, ",", y_coord, ")"),
    x = "Value",
    y = "Height (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  ) +
  scale_y_continuous(breaks = seq(0, 60, by = 20))

print(faceted_plot)

pdf("../../figs/mc_output/mc_plot_avg_metrics.pdf")
print(faceted_plot)
dev.off()