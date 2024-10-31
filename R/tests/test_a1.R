# Run this test from within the R/ directory to ensure
# correct interpretation of relative paths within the script.

dir_exp <- "tests/data/microhabitat"
dir_rec <- "tests/data/output_a1"

# Test dimPlot
dimPlot_exp <- readRDS(file.path(dir_exp, "dimPlot.rds"))
dimPlot_rec <- readRDS(file.path(dir_rec, "dimPlot.rds"))

stopifnot(identical(length(dimPlot_exp), length(dimPlot_rec)))
stopifnot(identical(dimPlot_exp, dimPlot_rec))

# Test Microhabitat
Nfiles <- 3
for (t in seq_len(Nfiles)) {
    microhab_filename <- paste("MicrohabitatMatrix", t, ".rds", sep="")

    microhab_exp <- readRDS(file.path(dir_exp, microhab_filename))
    microhab_rec <- readRDS(file.path(dir_rec, microhab_filename))

    stopifnot(identical(dim(microhab_exp), dim(microhab_rec)))
    stopifnot(isTRUE(all.equal(microhab_exp, microhab_rec, tolerance=1e-5)))
}
