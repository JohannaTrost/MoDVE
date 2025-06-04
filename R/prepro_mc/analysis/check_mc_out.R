library(microclimf)
library(terra)
library(readr)
library(viridis)
library(ggplot2)
library(lubridate)
library(reshape2)
library(dplyr)

# Function to load tair and average over time
load_avg_tair <- function(x, y, base_path, var = "tair") {
  file <- sprintf("%s/mc_x%d_y%d_v1.rds", base_path, x, y)
  data <- readRDS(file)
  tair_avg <- rowMeans(data[[var]], na.rm = FALSE)  # Average over time
  return(tair_avg)
}

# Settings
base_path <- "/Users/johanna/Uni/masterarbeit/data/mc_output"
z_levels <- 59
x_range <- 1:3
y_range <- 1:50

# Sample random slices
set.seed(42)
random_x <- sample(x_range, 3)  # 3 random X values for YZ slices
random_y <- sample(y_range, 3)  # 3 random Y values for XZ slices

# -------------- Check MC in XZ and YZ slices --------------

for (var in c("tair", "tcanopy", "relhum", "windspeed")) {
  # Gather data for XZ slices
  xz_data <- do.call(rbind, lapply(random_y, function(y) {
    do.call(rbind, lapply(x_range, function(x) {
      tair_avg <- load_avg_tair(x, y, base_path, var)
      data.frame(x = x, z = 1:z_levels, value = tair_avg, slice = paste("y =", y), variable = var)
    }))
  }))

  # Gather data for YZ slices
  yz_data <- do.call(rbind, lapply(random_x, function(x) {
    do.call(rbind, lapply(y_range, function(y) {
      tair_avg <- load_avg_tair(x, y, base_path, var)
      data.frame(y = y, z = 1:z_levels, value = tair_avg, slice = paste("x =", x), variable = var)
    }))
  }))

  # Plot XZ slices
  ggplot(xz_data, aes(x = x, y = z, fill = value)) +
    geom_tile() +
    scale_fill_viridis(option = "viridis", name = paste("Avg", var)) +
    facet_wrap(~ slice + variable, ncol = 1) +
    labs(title = paste("XZ Slices of Average", var), x = "X", y = "Height (Z)") +
    theme_minimal() +
    ggsave(paste0("xz_slices_", var, ".png"))

  # Plot YZ slices
  ggplot(yz_data, aes(x = y, y = z, fill = value)) +
    geom_tile() +
    scale_fill_viridis(option = "magma", name = paste("Avg", var)) +
    facet_wrap(~ slice + variable, ncol = 1) +
    labs(title = paste("YZ Slices of Average", var), x = "Y", y = "Height (Z)") +
    theme_minimal() +
    ggsave(paste0("yz_slices_", var, ".png"))
}

# -------------- Check time series of avg MC for this one year --------------

# Load all data and convet from hourly to mean daily values (50 x 50 x 59 x 365) for each variable
dates_2024 <- seq.Date(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")

x_range <- 1:3
y_range <- 1:50
z_levels <- 59
n_days <- 366
vars <- c("tair", "tcanopy", "relhum", "windspeed")

# Initialize storage
data_daily <- list()
for (v in vars) {
  data_daily[[v]] <- array(NA, dim = c(50, 50, z_levels, n_days))
}

# Helper to get daily means from hourly
hourly_to_daily <- function(mat) {
  daily <- matrix(NA, nrow = nrow(mat), ncol = 366)
  for (i in 1:366) {
    idx <- ((i - 1) * 24 + 1):(i * 24)
    if (length(idx) <= ncol(mat)) {
      daily[, i] <- rowMeans(mat[, idx], na.rm = TRUE)
    }
  }
  return(daily)
}

# Load all data
for (x in x_range) {
  for (y in y_range) {
    file <- sprintf("%s/mc_x%d_y%d_v1.rds", base_path, x, y)
    data <- readRDS(file)

    for (v in vars) {
      daily_mat <- hourly_to_daily(data[[v]])
      data_daily[[v]][x, y, , ] <- daily_mat
    }
  }
}

# Convert to time series format for ggplot
plot_data <- list()

for (v in vars) {
  # Average over grid (x, y)
  mean_grid <- apply(data_daily[[v]], c(3, 4), mean, na.rm = TRUE)  # [z, day]

  df <- data.frame(
    date = rep(dates_2024, times = 3),
    value = c(mean_grid[2, ], mean_grid[15, ], mean_grid[50, ]),
    height = factor(rep(c("z = 1.5m", "z = 14.5m", "z = 49.5"), each = n_days))
  )
  df$variable <- v
  plot_data[[v]] <- df
}

combined_df <- bind_rows(plot_data)

# Clean temp
combined_df <- combined_df %>%
  mutate(value = ifelse(variable %in% c("tair", "tcanopy") & value > 60, NA, value))

# Plot time series for each variable
p1 <- ggplot(combined_df, aes(x = date, y = value, color = height)) +
  geom_line(size = 0.7, alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  labs(title = "Daily Mean Time Series (Grid Avg.)",
       x = "Date (2024)", y = "Value", color = "Height (z)") +
  theme_minimal() +
  scale_color_manual(values = c("z = 1.5m" = "darkred",
                                "z = 14.5m" = "steelblue",
                                "z = 49.5" = "forestgreen"))
# Print plots to a pdf file
pdf("../../figs/mc_output/ts_24_mc_v1.pdf")
print(p1)
dev.off()

# --- Rel Hum profile

# Average over x, y, and time (dims 1, 2, 4)
relhum_profile <- apply(data_daily[["relhum"]], 3, mean, na.rm = TRUE)  # length 59
z_levels <- 1:59

profile_df <- data.frame(
  height = z_levels,
  relhum = relhum_profile
)

ggplot(profile_df, aes(x = relhum, y = height)) +
  geom_line(color = "dodgerblue", size = 1.2) +
  scale_y_reverse() +  # Optional: reverse to show height increasing upwards
  labs(
    title = "Average Relative Humidity Profile",
    x = "Relative Humidity (%)",
    y = "Height (z-level)"
  ) +
  theme_minimal()
