#!/usr/bin/env Rscript

# Run independent mzML/mass-accuracy CSV pairs concurrently. The launcher only
# examines files directly inside the input directory. Each pair is placed in an
# isolated working directory with temporary input symlinks before HDPairFinder
# runs, preventing its setwd() calls and fixed output names from colliding. The
# temporary links are removed after each job finishes.
#
# Usage:
#   Rscript run_parallel.R [input_directory] [workers] [output_directory]
#
# Defaults:
#   input_directory  <repository>/mzML
#   workers          up to 3 concurrent jobs and 90% of available CPU cores
#   output_directory <input_directory>/.hdpairfinder_runs
#
# Resource limits can be adjusted with these environment variables:
#   HDPAIRFINDER_CPU_FRACTION      CPU capacity ceiling, 0-1 (default: 0.90)
#   HDPAIRFINDER_MAX_WORKERS       hard worker cap (default: 3)
#   HDPAIRFINDER_RESERVED_CORES    additional CPU cores kept free (default: 0)
#   HDPAIRFINDER_THREADS_PER_JOB   BLAS/OpenMP threads per job (default: 1)
#   HDPAIRFINDER_NICE              Unix process niceness, 0-19 (default: 10)
#   HDPAIRFINDER_TELEMETRY_SECONDS sampling interval (default: 5)
#   HDPAIRFINDER_SEGFAULT_RETRIES  serial retries after exit 139 (default: 1)

launcher_directory <- function() {
        file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
        if (length(file_argument) == 0) {
                return(normalizePath(getwd(), mustWork = TRUE))
        }
        dirname(normalizePath(sub("^--file=", "", file_argument[1]), mustWork = TRUE))
}

is_inside <- function(path, directory) {
        path <- normalizePath(path, mustWork = FALSE)
        directory <- normalizePath(directory, mustWork = FALSE)
        identical(path, directory) || startsWith(path, paste0(directory, .Platform$file.sep))
}

read_positive_integer <- function(value, label) {
        parsed <- suppressWarnings(as.integer(value))
        if (length(parsed) != 1 || is.na(parsed) || parsed < 1 || as.character(parsed) != value) {
                stop(label, " must be a positive integer; received: ", value)
        }
        parsed
}

read_nonnegative_integer <- function(value, label, maximum = Inf) {
        parsed <- suppressWarnings(as.integer(value))
        if (length(parsed) != 1 || is.na(parsed) || parsed < 0 ||
            parsed > maximum || as.character(parsed) != value) {
                stop(
                        label, " must be an integer between 0 and ", maximum,
                        "; received: ", value
                )
        }
        parsed
}

read_fraction <- function(value, label) {
        parsed <- suppressWarnings(as.numeric(value))
        if (length(parsed) != 1 || !is.finite(parsed) || parsed <= 0 || parsed > 1) {
                stop(label, " must be greater than 0 and no greater than 1; received: ", value)
        }
        parsed
}

read_proc_cpu <- function() {
        if (!file.exists("/proc/stat")) {
                return(c(total = NA_real_, idle = NA_real_))
        }
        fields <- strsplit(trimws(readLines("/proc/stat", n = 1L, warn = FALSE)), "\\s+")[[1]]
        values <- suppressWarnings(as.numeric(fields[-1]))
        if (length(values) < 4 || anyNA(values)) {
                return(c(total = NA_real_, idle = NA_real_))
        }
        c(total = sum(values), idle = values[4] + if (length(values) >= 5) values[5] else 0)
}

read_proc_memory <- function() {
        if (!file.exists("/proc/meminfo")) {
                return(c(
                        total_gb = NA_real_, used_gb = NA_real_, available_gb = NA_real_,
                        used_percent = NA_real_, swap_used_gb = NA_real_
                ))
        }
        lines <- readLines("/proc/meminfo", warn = FALSE)
        value_kb <- function(key) {
                line <- lines[startsWith(lines, paste0(key, ":"))]
                if (length(line) == 0) return(NA_real_)
                suppressWarnings(as.numeric(strsplit(trimws(line[1]), "\\s+")[[1]][2]))
        }
        total <- value_kb("MemTotal")
        available <- value_kb("MemAvailable")
        swap_total <- value_kb("SwapTotal")
        swap_free <- value_kb("SwapFree")
        used <- total - available
        kb_per_gb <- 1024^2
        c(
                total_gb = total / kb_per_gb,
                used_gb = used / kb_per_gb,
                available_gb = available / kb_per_gb,
                used_percent = 100 * used / total,
                swap_used_gb = (swap_total - swap_free) / kb_per_gb
        )
}

sample_system_telemetry <- function(previous_cpu, started, active_jobs, detected_cores) {
        current_cpu <- read_proc_cpu()
        cpu_percent <- NA_real_
        if (!anyNA(c(previous_cpu, current_cpu))) {
                total_delta <- current_cpu["total"] - previous_cpu["total"]
                idle_delta <- current_cpu["idle"] - previous_cpu["idle"]
                if (total_delta > 0) cpu_percent <- 100 * (1 - idle_delta / total_delta)
        }

        load_average <- rep(NA_real_, 3)
        if (file.exists("/proc/loadavg")) {
                fields <- strsplit(readLines("/proc/loadavg", n = 1L, warn = FALSE), "\\s+")[[1]]
                load_average <- suppressWarnings(as.numeric(fields[1:3]))
        }
        memory <- read_proc_memory()
        row <- data.frame(
                timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
                active_jobs = active_jobs,
                system_cpu_percent = cpu_percent,
                load_1m = load_average[1],
                load_5m = load_average[2],
                load_15m = load_average[3],
                load_1m_percent_of_cores = 100 * load_average[1] / detected_cores,
                memory_total_gb = unname(memory["total_gb"]),
                memory_used_gb = unname(memory["used_gb"]),
                memory_available_gb = unname(memory["available_gb"]),
                memory_used_percent = unname(memory["used_percent"]),
                swap_used_gb = unname(memory["swap_used_gb"]),
                stringsAsFactors = FALSE
        )
        list(row = row, cpu = current_cpu)
}

append_telemetry <- function(row, telemetry_file) {
        write.table(
                row,
                telemetry_file,
                sep = ",",
                row.names = FALSE,
                col.names = !file.exists(telemetry_file),
                append = file.exists(telemetry_file),
                quote = TRUE
        )
}

safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)

discover_jobs <- function(input_directory, output_directory) {
        mzml_files <- list.files(
                input_directory,
                pattern = "\\.mzML$",
                recursive = FALSE,
                full.names = TRUE,
                ignore.case = TRUE
        )
        mzml_files <- mzml_files[!vapply(
                mzml_files,
                is_inside,
                logical(1),
                directory = output_directory
        )]

        csv_files <- sub("\\.mzML$", "_mass_accuracy.csv", mzml_files, ignore.case = TRUE)
        paired <- file.exists(csv_files)
        if (any(!paired)) {
                warning(
                        "Skipping ", sum(!paired), " mzML file(s) without a matching ",
                        "<mzML stem>_mass_accuracy.csv file:\n  ",
                        paste(mzml_files[!paired], collapse = "\n  "),
                        call. = FALSE
                )
        }

        mzml_files <- normalizePath(mzml_files[paired], mustWork = TRUE)
        csv_files <- normalizePath(csv_files[paired], mustWork = TRUE)
        if (length(mzml_files) == 0) {
                stop("No complete mzML/mass-accuracy CSV pairs were found in ", input_directory)
        }

        sample_ids <- tools::file_path_sans_ext(basename(mzml_files))
        sample_ids <- gsub("[^A-Za-z0-9._-]", "_", sample_ids)
        sample_ids <- make.unique(sample_ids, sep = "_")

        Map(
                function(sample_id, mzml, csv) {
                        list(sample_id = sample_id, mzml = mzml, csv = csv)
                },
                sample_ids,
                mzml_files,
                csv_files
        )
}

prepare_job <- function(job, run_directory, database_file = NULL) {
        job$working_directory <- file.path(run_directory, job$sample_id)
        if (!dir.create(job$working_directory, recursive = TRUE, showWarnings = FALSE)) {
                stop("Could not create job directory: ", job$working_directory)
        }

        link_sources <- c(job$mzml, job$csv)
        if (!is.null(database_file) && file.exists(database_file)) {
                link_sources <- c(link_sources, database_file)
        }
        link_targets <- file.path(job$working_directory, basename(link_sources))
        linked <- file.symlink(link_sources, link_targets)
        if (!all(linked)) {
                unlink(link_targets[linked])
                stop("Could not create input symlink(s) for job ", job$sample_id)
        }

        job$input_links <- link_targets
        job$log_file <- file.path(job$working_directory, "HDPairFinder.log")
        job$resource_file <- file.path(job$working_directory, "resources.txt")
        job
}

read_job_resources <- function(resource_file) {
        empty <- c(
                user_cpu_seconds = NA_real_, system_cpu_seconds = NA_real_,
                average_cpu_percent = NA_real_, max_rss_mb = NA_real_
        )
        if (!file.exists(resource_file)) return(empty)
        lines <- trimws(readLines(resource_file, warn = FALSE))
        extract_value <- function(label) {
                prefix <- paste0(label, ":")
                line <- lines[startsWith(lines, prefix)]
                if (length(line) == 0) return(NA_character_)
                trimws(substring(line[1], nchar(prefix) + 1L))
        }
        user_cpu <- suppressWarnings(as.numeric(extract_value("User time (seconds)")))
        system_cpu <- suppressWarnings(as.numeric(extract_value("System time (seconds)")))
        average_cpu <- suppressWarnings(as.numeric(sub("%$", "", extract_value("Percent of CPU this job got"))))
        max_rss_kb <- suppressWarnings(as.numeric(extract_value("Maximum resident set size (kbytes)")))
        c(
                user_cpu_seconds = user_cpu,
                system_cpu_seconds = system_cpu,
                average_cpu_percent = average_cpu,
                max_rss_mb = max_rss_kb / 1024
        )
}

run_job <- function(job, worker_script, threads_per_job, nice_value, attempt = 1L) {
        started <- Sys.time()
        log_file <- if (attempt == 1L) {
                job$log_file
        } else {
                file.path(
                        job$working_directory,
                        paste0("HDPairFinder_attempt_", attempt, ".log")
                )
        }
        resource_file <- if (attempt == 1L) {
                job$resource_file
        } else {
                file.path(
                        job$working_directory,
                        paste0("resources_attempt_", attempt, ".txt")
                )
        }
        rscript <- file.path(R.home("bin"), "Rscript")
        command <- rscript
        command_arguments <- c("--vanilla", shQuote(worker_script), shQuote(job$working_directory))

        nice_command <- Sys.which("nice")
        if (.Platform$OS.type == "unix" && nice_value > 0 && nzchar(nice_command)) {
                command <- nice_command
                command_arguments <- c(
                        "-n", as.character(nice_value), shQuote(rscript), command_arguments
                )
        }

        thread_environment <- paste0(
                c(
                        "OMP_NUM_THREADS=",
                        "OPENBLAS_NUM_THREADS=",
                        "MKL_NUM_THREADS=",
                        "VECLIB_MAXIMUM_THREADS=",
                        "NUMEXPR_NUM_THREADS="
                ),
                threads_per_job
        )
        thread_environment <- c(thread_environment, "LC_ALL=C")

        time_command <- Sys.which("time")
        if (.Platform$OS.type == "unix" && nzchar(time_command)) {
                timed_command <- command
                command <- time_command
                command_arguments <- c(
                        "-v", "-o", shQuote(resource_file),
                        shQuote(timed_command), command_arguments
                )
        }
        status <- suppressWarnings(system2(
                command = command,
                args = command_arguments,
                stdout = log_file,
                stderr = log_file,
                env = thread_environment
        ))
        resources <- read_job_resources(resource_file)

        list(
                sample_id = job$sample_id,
                attempt = attempt,
                success = identical(as.integer(status), 0L),
                status = as.integer(status),
                elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
                user_cpu_seconds = unname(resources["user_cpu_seconds"]),
                system_cpu_seconds = unname(resources["system_cpu_seconds"]),
                average_cpu_percent = unname(resources["average_cpu_percent"]),
                max_rss_mb = unname(resources["max_rss_mb"]),
                working_directory = job$working_directory,
                log_file = log_file,
                resource_file = resource_file
        )
}

main <- function() {
        if (!requireNamespace("future", quietly = TRUE)) {
                stop(
                        "The 'future' package is required. Install it with ",
                        "install.packages(\"future\")."
                )
        }

        launcher_dir <- launcher_directory()
        worker_script <- normalizePath(
                file.path(launcher_dir, "run_hdpairfinder_serial.R"),
                mustWork = TRUE
        )
        arguments <- commandArgs(trailingOnly = TRUE)
        input_directory <- if (length(arguments) >= 1) arguments[1] else file.path(launcher_dir, "mzML")
        input_directory <- normalizePath(input_directory, mustWork = TRUE)
        output_directory <- if (length(arguments) >= 3) {
                arguments[3]
        } else {
                file.path(input_directory, ".hdpairfinder_runs")
        }
        output_directory <- normalizePath(output_directory, mustWork = FALSE)

        jobs <- discover_jobs(input_directory, output_directory)
        detected_cores <- as.integer(parallelly::availableCores())
        cpu_fraction <- read_fraction(
                Sys.getenv("HDPAIRFINDER_CPU_FRACTION", unset = "0.90"),
                "HDPAIRFINDER_CPU_FRACTION"
        )
        reserved_cores <- read_nonnegative_integer(
                Sys.getenv("HDPAIRFINDER_RESERVED_CORES", unset = "0"),
                "HDPAIRFINDER_RESERVED_CORES"
        )
        threads_per_job <- read_positive_integer(
                Sys.getenv("HDPAIRFINDER_THREADS_PER_JOB", unset = "1"),
                "HDPAIRFINDER_THREADS_PER_JOB"
        )
        nice_value <- read_nonnegative_integer(
                Sys.getenv("HDPAIRFINDER_NICE", unset = "10"),
                "HDPAIRFINDER_NICE",
                maximum = 19L
        )
        telemetry_seconds <- read_positive_integer(
                Sys.getenv("HDPAIRFINDER_TELEMETRY_SECONDS", unset = "5"),
                "HDPAIRFINDER_TELEMETRY_SECONDS"
        )
        segfault_retries <- read_nonnegative_integer(
                Sys.getenv("HDPAIRFINDER_SEGFAULT_RETRIES", unset = "1"),
                "HDPAIRFINDER_SEGFAULT_RETRIES"
        )

        fraction_core_limit <- max(1L, floor(detected_cores * cpu_fraction))
        reserved_core_limit <- max(1L, detected_cores - reserved_cores)
        cpu_budget <- min(fraction_core_limit, reserved_core_limit)
        cpu_worker_limit <- max(1L, cpu_budget %/% threads_per_job)
        max_workers <- read_positive_integer(
                Sys.getenv("HDPAIRFINDER_MAX_WORKERS", unset = "3"),
                "HDPAIRFINDER_MAX_WORKERS"
        )
        worker_limit <- max(1L, min(length(jobs), max_workers, cpu_worker_limit))
        requested_workers <- if (length(arguments) >= 2) {
                read_positive_integer(arguments[2], "workers")
        } else {
                worker_limit
        }
        workers <- min(requested_workers, worker_limit)
        if (requested_workers > workers) {
                warning(
                        "Requested ", requested_workers, " workers, but resource limits allow ",
                        workers, ". Adjust HDPAIRFINDER_CPU_FRACTION, HDPAIRFINDER_MAX_WORKERS, ",
                        "or HDPAIRFINDER_RESERVED_CORES only if the system has enough CPU and memory.",
                        call. = FALSE
                )
        }

        if (!dir.exists(output_directory) &&
            !dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)) {
                stop("Could not create output directory: ", output_directory)
        }
        run_directory <- tempfile(pattern = "run_", tmpdir = output_directory)
        if (!dir.create(run_directory, showWarnings = FALSE)) {
                stop("Could not create run directory: ", run_directory)
        }

        database_file <- file.path(launcher_dir, "AMINES_library.csv")
        jobs <- lapply(jobs, prepare_job, run_directory = run_directory, database_file = database_file)
        on.exit(unlink(unlist(lapply(jobs, `[[`, "input_links"))), add = TRUE)
        manifest <- data.frame(
                sample_id = vapply(jobs, `[[`, character(1), "sample_id"),
                mzml = vapply(jobs, `[[`, character(1), "mzml"),
                mass_accuracy_csv = vapply(jobs, `[[`, character(1), "csv"),
                working_directory = vapply(jobs, `[[`, character(1), "working_directory"),
                stringsAsFactors = FALSE
        )
        write.csv(manifest, file.path(run_directory, "manifest.csv"), row.names = FALSE)

        initial_memory <- read_proc_memory()
        settings <- data.frame(
                detected_cores = detected_cores,
                cpu_fraction = cpu_fraction,
                cpu_core_budget = cpu_budget,
                threads_per_job = threads_per_job,
                inner_biocparallel = "SerialParam",
                requested_workers = requested_workers,
                workers = workers,
                hard_worker_cap = max_workers,
                reserved_cores = reserved_cores,
                nice = nice_value,
                telemetry_interval_seconds = telemetry_seconds,
                segfault_retries = segfault_retries,
                memory_total_gb = unname(initial_memory["total_gb"]),
                stringsAsFactors = FALSE
        )
        write.csv(settings, file.path(run_directory, "settings.csv"), row.names = FALSE)

        message("Processing ", length(jobs), " pair(s) with ", workers, " worker(s).")
        message(
                "Resource limits: ", threads_per_job, " thread(s)/job, ",
                round(cpu_fraction * 100), "% CPU ceiling (", cpu_budget, "/", detected_cores,
                " cores), worker cap ", max_workers,
                if (.Platform$OS.type == "unix") paste0(", nice ", nice_value) else ""
        )
        message("Run directory: ", run_directory)
        message("Inner xcms parallelism: BiocParallel SerialParam")
        message("Telemetry: ", file.path(run_directory, "telemetry.csv"))

        previous_plan <- future::plan()
        on.exit(future::plan(previous_plan), add = TRUE)
        future::plan(future::multisession, workers = workers)

        futures <- lapply(jobs, function(job) {
                future::future(
                        run_job(job, worker_script, threads_per_job, nice_value),
                        seed = TRUE
                )
        })
        telemetry_file <- file.path(run_directory, "telemetry.csv")
        telemetry_started <- Sys.time()
        previous_cpu <- c(total = NA_real_, idle = NA_real_)
        repeat {
                resolved <- vapply(futures, future::resolved, logical(1))
                telemetry <- sample_system_telemetry(
                        previous_cpu = previous_cpu,
                        started = telemetry_started,
                        active_jobs = sum(!resolved),
                        detected_cores = detected_cores
                )
                append_telemetry(telemetry$row, telemetry_file)
                previous_cpu <- telemetry$cpu
                if (all(resolved)) break
                Sys.sleep(telemetry_seconds)
        }
        results <- lapply(futures, future::value)

        if (segfault_retries > 0L) {
                for (retry_number in seq_len(segfault_retries)) {
                        retry_indexes <- which(vapply(
                                results,
                                function(result) identical(result$status, 139L),
                                logical(1)
                        ))
                        if (length(retry_indexes) == 0L) break

                        message(
                                "Retrying ", length(retry_indexes),
                                " job(s) that exited with SIGSEGV (status 139) serially; ",
                                "retry ", retry_number, "/", segfault_retries, "."
                        )
                        for (job_index in retry_indexes) {
                                results[[job_index]] <- run_job(
                                        jobs[[job_index]],
                                        worker_script,
                                        threads_per_job,
                                        nice_value,
                                        attempt = retry_number + 1L
                                )
                        }
                }
        }

        summary <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
        write.csv(summary, file.path(run_directory, "summary.csv"), row.names = FALSE)
        print(
                summary[, c(
                        "sample_id", "success", "status", "elapsed_seconds",
                        "attempt", "average_cpu_percent", "max_rss_mb", "log_file"
                )],
                row.names = FALSE
        )

        telemetry_data <- read.csv(telemetry_file)
        telemetry_summary <- data.frame(
                samples = nrow(telemetry_data),
                peak_system_cpu_percent = safe_max(telemetry_data$system_cpu_percent),
                peak_load_1m = safe_max(telemetry_data$load_1m),
                peak_load_1m_percent_of_cores = safe_max(telemetry_data$load_1m_percent_of_cores),
                peak_memory_used_gb = safe_max(telemetry_data$memory_used_gb),
                peak_memory_used_percent = safe_max(telemetry_data$memory_used_percent),
                minimum_memory_available_gb = safe_min(telemetry_data$memory_available_gb)
        )
        write.csv(
                telemetry_summary,
                file.path(run_directory, "telemetry_summary.csv"),
                row.names = FALSE
        )

        failures <- summary[!summary$success, , drop = FALSE]
        if (nrow(failures) > 0) {
                stop(
                        nrow(failures), " job(s) failed. Inspect the logs listed in ",
                        file.path(run_directory, "summary.csv"),
                        call. = FALSE
                )
        }

        message("All jobs completed successfully.")
}

main()
