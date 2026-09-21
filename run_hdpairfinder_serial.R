#!/usr/bin/env Rscript

# Run the upstream HDPairFinder script with serial BiocParallel execution.
# Independent samples are already parallelized by run_parallel.R, so allowing
# xcmsSet() to fork again adds nested workers and can make R unstable.

wrapper_directory <- function() {
        file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
        if (length(file_argument) == 0) {
                return(normalizePath(getwd(), mustWork = TRUE))
        }
        dirname(normalizePath(sub("^--file=", "", file_argument[1]), mustWork = TRUE))
}

register_serial_biocparallel <- function() {
        if (!requireNamespace("BiocParallel", quietly = TRUE)) {
                stop("The 'BiocParallel' package is required by HDPairFinder.")
        }
        parameter <- BiocParallel::SerialParam()
        BiocParallel::register(parameter, default = TRUE)
        invisible(parameter)
}

main <- function() {
        main_script <- normalizePath(
                file.path(wrapper_directory(), "HDPairFinder_v1.R"),
                mustWork = TRUE
        )
        register_serial_biocparallel()
        message("BiocParallel backend: SerialParam (inner xcms parallelism disabled).")
        source(main_script, chdir = FALSE)
}

if (sys.nframe() == 0L) main()
