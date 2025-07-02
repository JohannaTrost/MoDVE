# Prior to running this test, execute the A4 script using the input file
# located in this directory. This will create the necessary output files,
# which will subsequently be compared against the expected files.
# Run this test from within the R/ directory to ensure
# correct interpretation of relative paths within the script.

dir_expected <- "tests/data/growth_simulation"
dir_received <- "data/output_a5"

# Test simulation output files
timesteps <- c(2, 3)
spec_pool_id <- c(4, 5)
replicate <- 1
for (i in spec_pool_id) {
    sim_directory <- paste("ID_SpeciesP_", i, "_Rep_", replicate, sep="")

    # Test community summary file
    commun_expected <- read.csv(file.path(dir_expected, sim_directory, "CommunitySummary.csv"))
    commun_received <- read.csv(file.path(dir_received, sim_directory, "CommunitySummary.csv"))
    for (colname in colnames(commun_expected)) {
        stopifnot(isTRUE(all.equal(commun_expected[[colname]], commun_received[[colname]], tolerance=1e-5)))
    }

    # Test community summary file
    spec_summ_expected <- read.csv(file.path(dir_expected, sim_directory, "SpeciesSummary.csv"))
    spec_summ_received <- read.csv(file.path(dir_received, sim_directory, "SpeciesSummary.csv"))
    for (colname in colnames(spec_summ_expected)) {
        stopifnot(isTRUE(all.equal(spec_summ_expected[[colname]], spec_summ_received[[colname]], tolerance=1e-5)))
    }

    # Test individual matrices
    for (j in timesteps) {
        indiv_filename <- paste("IndividualMatrixTimeStep", j, ".csv", sep="")

        indiv_expected <- read.csv(file.path(dir_expected, sim_directory, indiv_filename))
        indiv_received <- read.csv(file.path(dir_received, sim_directory, indiv_filename))
        for (colname in colnames(indiv_expected)) {
            stopifnot(isTRUE(all.equal(indiv_expected[[colname]], indiv_received[[colname]], tolerance=1e-5)))
        }
    }
}
