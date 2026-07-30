# -----
# Useful functions for modeling pipeline

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


int_seq <- function(from, to, by=1L) {
    #' Custom sequence function for integers and floats with zero fractional part
    #'
    #' This function generates a sequence of integers similar to `seq.int`, but
    #' with additional checks and behaviors:
    #'
    #' * It avoids using the colon operator (':').
    #' * It ensures that `from`, `to` and `by` are integers or floats with zero fractional part.
    #' * It ensures that `by` is non-zero.
    #'
    #' @param from The starting value of the sequence (integer or float with zero fractional part).
    #' @param to The ending value of the sequence (integer or float with zero fractional part).
    #' @param by The increment between values in the sequence (non-zero integer or float with zero
    #'           fractional part, default=1).
    #'
    #' @return
    #' * A vector of integers generated using `seq.int` if the input is valid.
    #' * An empty vector if `by` is positive and `from` is greater than `to`,
    #'   or `by` is negative and `from` is less than `to`.
    #' * An error message if either `from`, `to` or `by` has a non-zero fractional part,
    #'   or if `by` is zero.
    #'
    #' @examples
    #' # Valid usage:
    #' int_seq(from = 2, to = 5, by = 1)     # [1] 2 3 4 5
    #' int_seq(from = 5.0, to = 3, by = -1)  # [1] 5 4 3
    #' int_seq(from = 4, to = 4, by = -1)    # [1] 4
    #' int_seq(from = 4, to = 5, by = -1)    # integer(0)
    #' int_seq(from = 2, to = 1, by = 2)     # integer(0)
    #'
    #' # Invalid usage:
    #' int_seq(from = 2.5, to = 5, by = 1)   # Error: Fractional part in 'from', 'to' or 'by'
    #' int_seq(from = 2, to = 5, by = 0)     # Error: 'by' must be non-zero

    # Stop if fractional part is non-zero
    if (!(from %% 1 == 0) || !(to %% 1 == 0) || !(by %% 1 == 0)) {
        stop("int_seq: Fractional part in 'from', 'to' or 'by'")
    }

    # Stop if `by` is zero
    if (by == 0) {
        stop("int_seq: 'by' must be non-zero")
    }

    # Convert to integers
    from <- as.integer(from)
    to <- as.integer(to)
    by <- as.integer(by)

    # Return empty integer array if `from` greater than `to` and `by` positive
    # or `from` less than `to` and `by` negative.
    if ((from > to && by > 0) || (from < to && by < 0)) {
        s <- integer()
    } else {
        s <- seq.int(from=from, to=to, by=by)
    }

    return(s)
}
