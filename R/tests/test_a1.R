# Prior to running this test, execute the A1 script using the input file
# located in this directory. This will create the necessary output files,
# which will subsequently be compared against the expected files.
# Run this test from within the R/ directory to ensure
# correct interpretation of relative paths within the script.

dir_expected <- "tests/data/microhabitat"
dir_received <- "tests/data/output_a1"

# Test dimPlot
dimPlot_expected <- readRDS(file.path(dir_expected, "dimPlot.rds"))
dimPlot_received <- readRDS(file.path(dir_received, "dimPlot.rds"))

stopifnot(identical(length(dimPlot_expected), length(dimPlot_received)))
stopifnot(identical(dimPlot_expected, dimPlot_received))

# Test microhabitat
Nfiles <- 3
for (t in seq_len(Nfiles)) {
    microhab_filename <- paste("microhabitatMatrix", t, ".rds", sep="")

    microhab_expected <- readRDS(file.path(dir_expected, microhab_filename))
    microhab_received <- readRDS(file.path(dir_received, microhab_filename))

    stopifnot(identical(dim(microhab_expected), dim(microhab_received)))
    stopifnot(isTRUE(all.equal(microhab_expected, microhab_received, tolerance=1e-5)))
}
