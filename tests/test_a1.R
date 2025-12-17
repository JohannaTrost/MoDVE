# Prior to running this test, execute the A1 script using the input file
# located in this directory. This will create the necessary output files,
# which will subsequently be compared against the expected files.
# Run this test from within the R/ directory to ensure
# correct interpretation of relative paths within the script.

dir_expected <- "tests/data/microhabitat"
dir_received <- "tests/data/output_a1"

# Test Microhabitat
Nfiles <- 3
t <- 3
for (t in seq_len(Nfiles)) {
    microhab_filename <- paste("MicrohabitatMatrix", t, ".rds", sep="")

    microhab_expected <- readRDS(file.path(dir_expected, microhab_filename))
    microhab_received <- readRDS(file.path(dir_received, paste("MicrohabitatMatrix_new", t, ".rds", sep="")))

    stopifnot(identical(dim(microhab_expected), dim(microhab_received)))
    stopifnot(isTRUE(all.equal(microhab_expected, microhab_received, tolerance=1e-5)))
}

# Convert microhabitat to data frame
microhab_to_df <- function(microhab_mat) {
  dims <- dim(microhab_mat)
  microhab_df <- tidyr::expand_grid(x = seq_len(dims[1]),
                                    y = seq_len(dims[2]),
                                    z = seq_len(dims[3]))

  sa <- sa_loss <- light <- vector("numeric", prod(dims[1:3]))

  i <- 1
  for (x in seq_len(dims[1])) {
    for (y in seq_len(dims[2])) {
      for (z in seq_len(dims[3])) {
        sa[i] <- microhab_mat[x, y, z, 1]
        sa_loss[i] <- microhab_mat[x, y, z, 2]
        light[i] <- microhab_mat[x, y, z, 3]
        i <- i + 1
      }
    }
  }
  microhab_df$sa <- sa
  microhab_df$sa_loss <- sa_loss
  microhab_df$light <- light
  return(microhab_df)
}

df_expected <- microhab_to_df(microhab_expected)
df_received <- microhab_to_df(microhab_received)

df_expected$sa_loss[is.nan(df_expected$sa_loss)] <- 0
sum(df_expected$sa)
sum(df_received$sa)
sum(df_expected$sa) / sum(df_received$sa)
df_expected$sa - df_received$sa

df_diff <- tibble::tibble(
  "x" = df_expected$x,
  "y" = df_expected$y,
  "z" = df_expected$z,
  "exptd" = df_expected$sa,
  "recd" = df_received$sa,
  "diff" = exptd - recd
) |>
  dplyr::filter(diff != 0)

# don't drop voxels in corridor?

# Light is correct
any(abs(df_expected$light - df_received$light) > 0.000001)

sum(microhab_mat[,,,1])
sum(Mat_surface_per_cell[21:50, 21:50,])

sum(microhab_expected[,,,1])
sum(microhab_received[,,,1])
