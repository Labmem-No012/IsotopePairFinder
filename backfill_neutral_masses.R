#!/usr/bin/env Rscript

# Add unlabeled neutral [M] masses to existing HDPairFinder V1 CSV outputs.
#
# Usage:
#   Rscript backfill_neutral_masses.R INPUT_PATH [OUTPUT_DIRECTORY] [--overwrite]
#
# INPUT_PATH can be one CSV file or a directory. Directories are searched
# recursively. The original files are never modified; converted copies retain
# their relative paths under OUTPUT_DIRECTORY. When OUTPUT_DIRECTORY is omitted,
# a directory named "neutral_mass_outputs" is created beside/inside INPUT_PATH.

mass_n14 <- 14.00307400443
mass_n15 <- 15.00010889888
n15_shift <- mass_n15 - mass_n14
proton_mass <- 1.007276

required_pair_columns <- c("mz_light", "mz_heavy", "num_tag")
neutral_mass_columns <- c("unlab_light_mz", "unlab_heavy_mz")

# Reuse the 14N/15N monoisotopic masses and isotope shift from HDPairFinder's
# pair-picking module, plus its proton mass. Despite the legacy *_mz column
# names, the calculated values are neutral monoisotopic [M] masses.
add_unlabeled_neutral_masses <- function(feature_table) {
        missing_columns <- setdiff(required_pair_columns, colnames(feature_table))
        if (length(missing_columns) != 0) {
                stop(
                        "Cannot calculate neutral masses; missing column(s): ",
                        paste(missing_columns, collapse = ", ")
                )
        }

        feature_table$unlab_light_mz <- feature_table$mz_light - proton_mass
        feature_table$unlab_heavy_mz <- feature_table$mz_heavy -
                feature_table$num_tag * n15_shift - proton_mass
        feature_table
}

usage <- function() {
        paste(
                "Usage:",
                "  Rscript backfill_neutral_masses.R INPUT_PATH [OUTPUT_DIRECTORY] [--overwrite]",
                "",
                "INPUT_PATH may be a CSV file or a directory searched recursively.",
                "Original V1 outputs are never modified.",
                sep = "\n"
        )
}

normalize_output_path <- function(path) {
        normalizePath(path, winslash = "/", mustWork = FALSE)
}

is_path_within <- function(path, directory) {
        normalized_path <- normalize_output_path(path)
        normalized_directory <- sub("/+$", "", normalize_output_path(directory))
        identical(normalized_path, normalized_directory) ||
                startsWith(normalized_path, paste0(normalized_directory, "/"))
}

validate_numeric_column <- function(values, column_name, input_file) {
        character_values <- trimws(as.character(values))
        numeric_values <- suppressWarnings(as.numeric(character_values))
        invalid <- !is.na(values) & nzchar(character_values) & is.na(numeric_values)
        if (any(invalid)) {
                examples <- unique(character_values[invalid])
                stop(
                        "Column ", column_name, " contains non-numeric value(s) in ",
                        input_file, ": ", paste(utils::head(examples, 3), collapse = ", ")
                )
        }
        numeric_values
}

read_csv_header <- function(input_file) {
        colnames(utils::read.csv(
                input_file,
                nrows = 0,
                check.names = FALSE,
                stringsAsFactors = FALSE
        ))
}

write_csv_atomically <- function(feature_table, output_file, overwrite) {
        output_parent <- dirname(output_file)
        if (!dir.exists(output_parent) &&
            !dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)) {
                stop("Could not create output directory: ", output_parent)
        }

        if (file.exists(output_file) && !overwrite) {
                return(FALSE)
        }

        temporary_file <- tempfile(
                pattern = paste0(".", basename(output_file), "."),
                tmpdir = output_parent
        )
        on.exit(unlink(temporary_file), add = TRUE)
        utils::write.csv(feature_table, temporary_file, row.names = FALSE, na = "NA")

        if (file.exists(output_file) && !file.remove(output_file)) {
                stop("Could not replace existing output: ", output_file)
        }
        if (!file.rename(temporary_file, output_file)) {
                stop("Could not move completed output into place: ", output_file)
        }
        TRUE
}

discover_csv_files <- function(input_path, output_directory) {
        if (!dir.exists(input_path)) {
                return(input_path)
        }

        input_files <- list.files(
                input_path,
                pattern = "\\.csv$",
                recursive = TRUE,
                full.names = TRUE,
                ignore.case = TRUE
        )
        if (!is_path_within(output_directory, input_path)) {
                return(input_files)
        }
        input_files[!vapply(
                input_files,
                is_path_within,
                logical(1),
                directory = output_directory
        )]
}

relative_input_path <- function(input_file, input_path) {
        if (!dir.exists(input_path)) {
                return(basename(input_file))
        }

        input_root <- paste0(sub("/+$", "", input_path), "/")
        substring(input_file, nchar(input_root) + 1L)
}

process_csv_file <- function(input_file, output_file, overwrite) {
        header <- read_csv_header(input_file)
        if (!all(required_pair_columns %in% header)) {
                return(list(
                        status = "skipped_not_pair_output",
                        rows = NA_integer_,
                        details = paste(
                                "Missing",
                                paste(setdiff(required_pair_columns, header), collapse = ", ")
                        )
                ))
        }

        if (file.exists(output_file) && !overwrite) {
                return(list(
                        status = "skipped_existing_output",
                        rows = NA_integer_,
                        details = "Use --overwrite to replace the generated copy"
                ))
        }

        feature_table <- utils::read.csv(
                input_file,
                check.names = FALSE,
                stringsAsFactors = FALSE
        )
        for (column_name in required_pair_columns) {
                feature_table[[column_name]] <- validate_numeric_column(
                        feature_table[[column_name]],
                        column_name,
                        input_file
                )
        }

        already_had_neutral_masses <- all(neutral_mass_columns %in% colnames(feature_table))
        feature_table <- add_unlabeled_neutral_masses(feature_table)
        write_csv_atomically(feature_table, output_file, overwrite = overwrite)

        list(
                status = if (already_had_neutral_masses) "recalculated" else "converted",
                rows = nrow(feature_table),
                details = ""
        )
}

backfill_neutral_masses <- function(input_path, output_directory = NULL, overwrite = FALSE) {
        input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)
        if (!dir.exists(input_path) && !grepl("\\.csv$", input_path, ignore.case = TRUE)) {
                stop("INPUT_PATH must be a CSV file or directory: ", input_path)
        }

        if (is.null(output_directory)) {
                output_parent <- if (dir.exists(input_path)) input_path else dirname(input_path)
                output_directory <- file.path(output_parent, "neutral_mass_outputs")
        }
        output_directory <- normalize_output_path(output_directory)

        if (dir.exists(input_path) && identical(input_path, output_directory)) {
                stop("OUTPUT_DIRECTORY must differ from INPUT_PATH to protect V1 outputs.")
        }

        input_files <- discover_csv_files(input_path, output_directory)
        if (length(input_files) == 0) {
                stop("No CSV files found under INPUT_PATH: ", input_path)
        }
        input_files <- sort(normalizePath(input_files, winslash = "/", mustWork = TRUE))

        if (!dir.exists(output_directory) &&
            !dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)) {
                stop("Could not create OUTPUT_DIRECTORY: ", output_directory)
        }

        results <- vector("list", length(input_files))
        for (index in seq_along(input_files)) {
                input_file <- input_files[index]
                relative_path <- relative_input_path(input_file, input_path)
                output_file <- file.path(output_directory, relative_path)

                result <- tryCatch(
                        process_csv_file(input_file, output_file, overwrite = overwrite),
                        error = function(condition) {
                                list(
                                        status = "error",
                                        rows = NA_integer_,
                                        details = conditionMessage(condition)
                                )
                        }
                )
                results[[index]] <- data.frame(
                        input_file = input_file,
                        output_file = if (result$status == "skipped_not_pair_output") "" else output_file,
                        status = result$status,
                        rows = result$rows,
                        details = result$details,
                        stringsAsFactors = FALSE
                )

                if (result$status %in% c("converted", "recalculated")) {
                        message("[", result$status, "] ", relative_path, " (", result$rows, " rows)")
                } else if (identical(result$status, "error")) {
                        message("[error] ", relative_path, ": ", result$details)
                }
        }

        summary_table <- do.call(rbind, results)
        summary_file <- file.path(output_directory, "neutral_mass_conversion_summary.csv")
        utils::write.csv(summary_table, summary_file, row.names = FALSE, na = "")

        status_counts <- table(summary_table$status)
        message("Summary: ", summary_file)
        message(paste(
                paste(names(status_counts), as.integer(status_counts), sep = "="),
                collapse = "; "
        ))

        attr(summary_table, "summary_file") <- summary_file
        summary_table
}

main <- function() {
        arguments <- commandArgs(trailingOnly = TRUE)
        if (length(arguments) == 0 || any(arguments %in% c("-h", "--help"))) {
                cat(usage(), "\n")
                return(invisible(NULL))
        }

        overwrite <- "--overwrite" %in% arguments
        positional_arguments <- arguments[arguments != "--overwrite"]
        unknown_options <- positional_arguments[startsWith(positional_arguments, "-")]
        if (length(unknown_options) != 0) {
                stop("Unknown option(s): ", paste(unknown_options, collapse = ", "), "\n", usage())
        }
        if (length(positional_arguments) < 1 || length(positional_arguments) > 2) {
                stop(usage())
        }

        summary_table <- backfill_neutral_masses(
                input_path = positional_arguments[1],
                output_directory = if (length(positional_arguments) == 2) positional_arguments[2] else NULL,
                overwrite = overwrite
        )
        if (any(summary_table$status == "error")) {
                quit(save = "no", status = 1L)
        }
        invisible(summary_table)
}

if (sys.nframe() == 0L) main()
