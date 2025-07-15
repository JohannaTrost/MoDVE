

x_dim <- 50
y_dim <- 50
out_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_output"


# Initialize the microclimate matrix
max_hgt <- 59
n_temp_metrics <- 14
mc_matrix <- array(rep(NA, x_dim * y_dim * max_hgt * n_temp_metrics),
                   dim = c(x_dim, y_dim, max_hgt, n_temp_metrics))

# Populate the main matrix with results
total_time <- 0
successful_cells <- 0

for (x in 1:x_dim) {
  for (y in 1:y_dim) {
    # Create a unique filename for each cell
    file_path <- file.path(out_dir, paste0("/v3_mc_x", x, "_y", y, ".rds"))

    # Check if the file exists
    if (file.exists(file_path)) {
      # Load the data from the file
      result <- readRDS(file_path)

      # Ensure the result is not NULL
      if (!is.null(result)) {
        result$x <- x
        result$y <- y
        actual_heights <- nrow(result$data)

        # Fill the matrix (up to the actual number of heights)
        mc_matrix[x, y, 1:actual_heights, ] <- result$data

        numNAs <- sum(is.na(result$data))
        if (numNAs > 0) {
          cat("Warning: Found", numNAs, "NA values in cell (", x, ",", y, ").\n")
        }

        total_time <- total_time + result$processing_time
        successful_cells <- successful_cells + 1

        if (successful_cells %% 100 == 0) {
          avg_time <- total_time / successful_cells
          cat("Processed", successful_cells, "cells. Average time per cell:",
              round(avg_time, 3), "seconds\n")
        }
      }
    } else {
      stop(paste0("File not found: ", file_path))
    }
  }
}

cat("Processing complete! Total successful cells:", successful_cells, "\n")
cat("Total processing time:", round(total_time, 2), "seconds\n")
cat("Average time per cell:", round(total_time / successful_cells, 3), "seconds\n")

# Save the result
saveRDS(mc_matrix, "/Users/johanna/Uni/masterarbeit/data/mc_output/v4_2024_regua_mc_matrix.rds")

mc_test <- readRDS("/Users/johanna/Uni/masterarbeit/data/mc_output/v4_2024_regua_mc_matrix.rds")
