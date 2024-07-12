library("optparse")
library("configr")


parse_config <- function() {
    # Parse command line arguments
    parser <- OptionParser()
    parser <- add_option(parser,
        c("-i", "--input"),
        type="character",
        default=NA,
        metavar="PATH_TO_TOML",
        help="Path to TOML formatted input file")
    opt <- parse_args(parser)

    # Check that input file is provided
    if (is.na(opt$input)) {
        stop("Input parameter must be provided. See script usage (--help)")
    } else {
        filepath <- file.path(opt$input)
    }

    # Check whether provided input file exists and is in TOML format
    if (!(file.exists(filepath) && !dir.exists(filepath))) {
        stop("Input file doesn't exist")
    } else if (!is.toml.file(file=opt$input)) {
        stop("Input file isn't in TOML format")
    } else {
        config <- read.config(file=filepath)
    }

    return(config)
}
