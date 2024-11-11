# Prior to running this test, execute the A3 script using the input file
# located in this directory. This will create the necessary output files,
# which will subsequently be compared against the expected files.
# Run this test from within the R/ directory to ensure
# correct interpretation of relative paths within the script.

dir_expected <- "tests/data/epiphyte_distributions"
dir_received <- "tests/data/output_a3"

# Test distributed species pools
spec_pool_id <- c(4, 5)
replicate <- 1
for (i in spec_pool_id) {
    dist_filename <- paste("ID_SpeciesP_", i, "_Rep_", replicate, ".csv", sep="")

    dist_expected <- read.csv(file.path(dir_expected, dist_filename))
    dist_received <- read.csv(file.path(dir_received, dist_filename))

    # Test each column separately to allow reordering
    for (colname in colnames(dist_expected)) {
        stopifnot(isTRUE(all.equal(dist_expected[[colname]], dist_received[[colname]], tolerance=1e-5)))
    }
}
