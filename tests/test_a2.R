# Prior to running this test, execute the A2 script using the input file
# located in this directory. This will create the necessary output files,
# which will subsequently be compared against the expected files.
# Run this test from within the R/ directory to ensure
# correct interpretation of relative paths within the script.

dir_expected <- "tests/data/species_pools"
dir_received <- "tests/data/output_a2"

# Test Trait ranges
trait_expected <- read.table(file.path(dir_expected, "TraitRanges.csv"), sep=",", header=FALSE)
trait_received <- read.table(file.path(dir_received, "TraitRanges.csv"), sep=",", header=FALSE)

stopifnot(isTRUE(all.equal(trait_expected, trait_received, tolerance=1e-5)))

# Test Species pools
Nfiles <- 5
for (t in seq_len(Nfiles)) {
    spec_pool_filename <- paste("SpeciesPool", t, ".csv", sep="")

    spec_pool_expected <- read.csv(file.path(dir_expected, spec_pool_filename))
    spec_pool_received <- read.csv(file.path(dir_received, spec_pool_filename))

    # Test each column separately to allow reordering
    for (colname in colnames(spec_pool_expected)) {
        stopifnot(isTRUE(all.equal(spec_pool_expected[[colname]], spec_pool_received[[colname]], tolerance=1e-5)))
    }
}
