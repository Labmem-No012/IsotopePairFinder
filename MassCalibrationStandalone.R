# Standalone mass calibration for mzML files
#
# Edit the settings in this section, then run:
#   Rscript MassCalibrationStandalone.R

# User settings

# A single mzML file or a directory containing mzML files.
input_path <- "/home/haotian/IsotopePairFinder/mzML"

# Use NULL to write beside each input file, or provide another directory.
output_directory <- "/home/haotian/IsotopePairFinder/mzML/Mass_Accuracy"

positive_standards <- data.frame(
  theoretical_mz = c(145.0799, 161.0748, 258.0845, 260.0847), #167.0833,
  expected_rt_min = c(1.63038333333, 1.31658333333, 7.12723333333, 6.56471666667)# rep(NA_real_, 5L)
)

negative_standards <- data.frame(
  theoretical_mz = c(144.0398, 155.0652, 165.0687, 199.0329),
  expected_rt_min = rep(NA_real_, 4L)
)

eic_search_ppm <- 20
retention_time_tolerance_min <- 5
minimum_eic_intensity <- 1000
minimum_eic_scans <- 3L
minimum_eic_snr <- 5
expected_peak_width_seconds <- c(10, 120)

# The observed m/z is the intensity-weighted mean of the apex scan and one
# immediately adjacent MS1 scan on each side.
neighboring_scans_each_side <- 1L

generate_diagnostic_pdf <- TRUE
spectrum_plot_ppm <- 100

# Validation

check_required_packages <- function() {
  required_packages <- c("BiocParallel", "MSnbase", "mzR", "xcms")
  installed <- vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  missing_packages <- required_packages[!installed]

  if (length(missing_packages) > 0L) {
    install_expression <- paste(
      sprintf('"%s"', missing_packages),
      collapse = ", "
    )
    stop(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them with BiocManager::install(c(",
      install_expression,
      ")).",
      call. = FALSE
    )
  }
}

check_positive_number <- function(value, name, allow_zero = FALSE) {
  lower_bound_is_valid <- if (allow_zero) value >= 0 else value > 0
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        !lower_bound_is_valid) {
    comparison <- if (allow_zero) "non-negative" else "positive"
    stop(name, " must be one finite ", comparison, " number.", call. = FALSE)
  }
}

validate_standard_table <- function(standards, name) {
  required_columns <- c("theoretical_mz", "expected_rt_min")
  if (!is.data.frame(standards) ||
        !all(required_columns %in% names(standards))) {
    stop(
      name,
      " must be a data frame containing theoretical_mz and expected_rt_min.",
      call. = FALSE
    )
  }
  if (nrow(standards) == 0L) {
    stop(name, " must contain at least one standard.", call. = FALSE)
  }
  if (!is.numeric(standards$theoretical_mz) ||
        any(!is.finite(standards$theoretical_mz)) ||
        any(standards$theoretical_mz <= 0)) {
    stop(name, "$theoretical_mz must contain positive numbers.", call. = FALSE)
  }
  invalid_rt <- !is.na(standards$expected_rt_min) &
    (!is.finite(standards$expected_rt_min) | standards$expected_rt_min < 0)
  if (!is.numeric(standards$expected_rt_min) || any(invalid_rt)) {
    stop(
      name,
      "$expected_rt_min must contain non-negative minutes or NA.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

validate_configuration <- function(
    input_path,
    output_directory,
    positive_standards,
    negative_standards,
    settings) {
  if (length(input_path) != 1L || is.na(input_path) || !nzchar(input_path) ||
        !file.exists(input_path)) {
    stop(
      "input_path must be an existing mzML file or directory.",
      call. = FALSE
    )
  }
  invalid_output_directory <- length(output_directory) != 1L ||
    is.na(output_directory) || !nzchar(output_directory)
  if (!is.null(output_directory) && invalid_output_directory) {
    stop("output_directory must be NULL or one directory path.", call. = FALSE)
  }

  validate_standard_table(positive_standards, "positive_standards")
  validate_standard_table(negative_standards, "negative_standards")
  check_positive_number(settings$eic_search_ppm, "eic_search_ppm")
  check_positive_number(
    settings$retention_time_tolerance_min,
    "retention_time_tolerance_min"
  )
  check_positive_number(
    settings$minimum_eic_intensity,
    "minimum_eic_intensity",
    allow_zero = TRUE
  )
  check_positive_number(
    settings$minimum_eic_snr,
    "minimum_eic_snr",
    allow_zero = TRUE
  )
  check_positive_number(settings$spectrum_plot_ppm, "spectrum_plot_ppm")

  if (length(settings$minimum_eic_scans) != 1L ||
        !is.numeric(settings$minimum_eic_scans) ||
        !is.finite(settings$minimum_eic_scans) ||
        settings$minimum_eic_scans < 1L ||
        settings$minimum_eic_scans != as.integer(settings$minimum_eic_scans)) {
    stop("minimum_eic_scans must be one positive integer.", call. = FALSE)
  }
  if (length(settings$expected_peak_width_seconds) != 2L ||
        !is.numeric(settings$expected_peak_width_seconds) ||
        any(!is.finite(settings$expected_peak_width_seconds)) ||
        any(settings$expected_peak_width_seconds <= 0) ||
        settings$expected_peak_width_seconds[1L] >=
          settings$expected_peak_width_seconds[2L]) {
    stop(
      paste(
        "expected_peak_width_seconds must contain two increasing",
        "positive numbers."
      ),
      call. = FALSE
    )
  }
  if (length(settings$neighboring_scans_each_side) != 1L ||
        !is.numeric(settings$neighboring_scans_each_side) ||
        !is.finite(settings$neighboring_scans_each_side) ||
        settings$neighboring_scans_each_side < 0L ||
        settings$neighboring_scans_each_side !=
          as.integer(settings$neighboring_scans_each_side)) {
    stop(
      "neighboring_scans_each_side must be one non-negative integer.",
      call. = FALSE
    )
  }
  if (length(settings$generate_diagnostic_pdf) != 1L ||
        is.na(settings$generate_diagnostic_pdf) ||
        !is.logical(settings$generate_diagnostic_pdf)) {
    stop("generate_diagnostic_pdf must be TRUE or FALSE.", call. = FALSE)
  }

  invisible(TRUE)
}

# Input helpers

find_mzml_files <- function(path) {
  if (!dir.exists(path)) {
    if (!grepl("[.]mzML$", path, ignore.case = TRUE)) {
      stop("The input file must have an .mzML extension.", call. = FALSE)
    }
    return(normalizePath(path, mustWork = TRUE))
  }

  files <- list.files(
    path,
    pattern = "[.]mzML$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(files) == 0L) {
    stop("No mzML files were found in: ", path, call. = FALSE)
  }

  sort(normalizePath(files, mustWork = TRUE))
}

result_directory_for_file <- function(file, output_directory) {
  directory <- if (is.null(output_directory)) {
    dirname(file)
  } else {
    output_directory
  }
  if (!dir.exists(directory) &&
        !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory: ", directory, call. = FALSE)
  }
  normalizePath(directory, mustWork = TRUE)
}

detect_polarity <- function(metadata, file) {
  polarities <- metadata$polarity[
    metadata$msLevel == 1L &
      !is.na(metadata$polarity) &
      metadata$polarity != 0
  ]
  polarities <- unique(sign(polarities))

  if (length(polarities) == 0L) {
    stop(
      "No MS1 polarity metadata were found in ",
      basename(file),
      ".",
      call. = FALSE
    )
  }
  if (length(polarities) > 1L) {
    stop(
      "Both positive and negative MS1 spectra were found in ",
      basename(file),
      ". Split mixed-polarity data before calibration.",
      call. = FALSE
    )
  }

  if (polarities == 1L) "positive" else "negative"
}

retention_time_range <- function(expected_rt_min, tolerance_min,
                                 full_range_seconds) {
  if (is.na(expected_rt_min)) {
    return(full_range_seconds)
  }

  requested <- c(
    expected_rt_min - tolerance_min,
    expected_rt_min + tolerance_min
  ) * 60
  overlap <- c(
    max(requested[1L], full_range_seconds[1L]),
    min(requested[2L], full_range_seconds[2L])
  )
  if (overlap[1L] > overlap[2L]) NULL else overlap
}

# Peak measurement

estimate_eic_snr <- function(intensity, retention_time_seconds, apex_index,
                             peak_width_seconds, minimum_scans,
                             minimum_intensity) {
  detected_peaks <- tryCatch(
    suppressWarnings(xcms::peaksWithCentWave(
      int = intensity,
      rt = retention_time_seconds,
      peakwidth = peak_width_seconds,
      snthresh = 0,
      prefilter = c(minimum_scans, minimum_intensity),
      noise = 0,
      firstBaselineCheck = TRUE
    )),
    error = function(error) NULL
  )
  if (is.null(detected_peaks) || nrow(detected_peaks) == 0L ||
        !"sn" %in% colnames(detected_peaks)) {
    return(NA_real_)
  }

  apex_time <- retention_time_seconds[apex_index]
  matching_rows <- which(
    detected_peaks[, "rtmin"] <= apex_time &
      detected_peaks[, "rtmax"] >= apex_time
  )
  if (length(matching_rows) == 0L) {
    matching_rows <- which.min(abs(detected_peaks[, "rt"] - apex_time))
  } else if (length(matching_rows) > 1L) {
    matching_rows <- matching_rows[
      which.min(abs(detected_peaks[matching_rows, "rt"] - apex_time))
    ]
  }

  as.numeric(detected_peaks[matching_rows[1L], "sn"])
}

collect_mass_points <- function(ms_connection, metadata, ms1_rows,
                                apex_rt_seconds, mz_range,
                                neighboring_scans_each_side,
                                minimum_intensity) {
  apex_position <- which.min(
    abs(metadata$retentionTime[ms1_rows] - apex_rt_seconds)
  )
  offsets <- seq.int(
    -neighboring_scans_each_side,
    neighboring_scans_each_side
  )
  positions <- apex_position + offsets
  if (any(positions < 1L | positions > length(ms1_rows))) {
    return(list(points = data.frame(), apex_spectrum = NULL))
  }

  points <- vector("list", length(positions))
  apex_spectrum <- NULL
  for (index in seq_along(positions)) {
    scan_row <- ms1_rows[positions[index]]
    spectrum <- mzR::peaks(ms_connection, scan_row)
    if (offsets[index] == 0L) apex_spectrum <- spectrum

    candidates <- which(
      is.finite(spectrum[, "mz"]) &
        is.finite(spectrum[, "intensity"]) &
        spectrum[, "mz"] >= mz_range[1L] &
        spectrum[, "mz"] <= mz_range[2L] &
        spectrum[, "intensity"] >= minimum_intensity
    )
    if (length(candidates) == 0L) next

    strongest <- candidates[
      which.max(spectrum[candidates, "intensity"])
    ]
    points[[index]] <- data.frame(
      scan_row = scan_row,
      retention_time_min = metadata$retentionTime[scan_row] / 60,
      mz = spectrum[strongest, "mz"],
      intensity = spectrum[strongest, "intensity"]
    )
  }

  points <- points[!vapply(points, is.null, logical(1))]
  points <- if (length(points) == 0L) data.frame() else do.call(rbind, points)
  list(points = points, apex_spectrum = apex_spectrum)
}

absolute_ppm_error <- function(theoretical_mz, observed_mz) {
  abs(observed_mz - theoretical_mz) / theoretical_mz * 1e6
}

empty_standard_result <- function(file, polarity, theoretical_mz,
                                  expected_rt_min) {
  data.frame(
    file = basename(file),
    polarity = polarity,
    theoretical_mz = theoretical_mz,
    expected_rt_min = expected_rt_min,
    observed_mz = NA_real_,
    absolute_ppm_difference = NA_real_,
    apex_retention_time_min = NA_real_,
    eic_snr = NA_real_,
    accepted = FALSE,
    status = "not measured",
    file_average_absolute_ppm_difference = NA_real_,
    stringsAsFactors = FALSE
  )
}

measure_standard <- function(file, standard, polarity, ms_data,
                             ms_connection, metadata, ms1_rows,
                             full_rt_range, settings) {
  theoretical_mz <- standard$theoretical_mz
  expected_rt_min <- standard$expected_rt_min
  result <- empty_standard_result(
    file,
    polarity,
    theoretical_mz,
    expected_rt_min
  )
  diagnostic <- list(
    eic = data.frame(),
    mass_points = data.frame(),
    apex_spectrum = NULL
  )

  rt_range <- retention_time_range(
    expected_rt_min,
    settings$retention_time_tolerance_min,
    full_rt_range
  )
  if (is.null(rt_range)) {
    result$status <- "expected RT is outside the acquisition range"
    return(list(result = result, diagnostic = diagnostic))
  }

  mz_tolerance <- theoretical_mz * settings$eic_search_ppm * 1e-6
  mz_range <- c(
    theoretical_mz - mz_tolerance,
    theoretical_mz + mz_tolerance
  )
  chromatogram <- MSnbase::chromatogram(
    ms_data,
    mz = matrix(mz_range, nrow = 1L),
    rt = matrix(rt_range, nrow = 1L),
    aggregationFun = "max",
    missing = 0,
    msLevel = 1L,
    BPPARAM = BiocParallel::SerialParam()
  )[1L, 1L]

  eic_intensity <- as.numeric(MSnbase::intensity(chromatogram))
  eic_rt_seconds <- as.numeric(MSnbase::rtime(chromatogram))
  valid_rows <- is.finite(eic_rt_seconds)
  eic_rt_seconds <- eic_rt_seconds[valid_rows]
  eic_intensity <- eic_intensity[valid_rows]
  eic_intensity[!is.finite(eic_intensity)] <- 0
  diagnostic$eic <- data.frame(
    retention_time_min = eic_rt_seconds / 60,
    intensity = eic_intensity,
    used_for_mass_average = FALSE
  )

  detected <- eic_intensity >= settings$minimum_eic_intensity
  if (sum(detected) < settings$minimum_eic_scans) {
    result$status <- "too few EIC scans above the intensity threshold"
    return(list(result = result, diagnostic = diagnostic))
  }

  apex_index <- which.max(ifelse(detected, eic_intensity, -Inf))
  apex_rt_seconds <- eic_rt_seconds[apex_index]
  result$apex_retention_time_min <- apex_rt_seconds / 60
  result$eic_snr <- estimate_eic_snr(
    eic_intensity,
    eic_rt_seconds,
    apex_index,
    settings$expected_peak_width_seconds,
    settings$minimum_eic_scans,
    settings$minimum_eic_intensity
  )

  if (is.na(result$eic_snr)) {
    result$status <- "EIC SNR could not be estimated"
    return(list(result = result, diagnostic = diagnostic))
  }
  if (result$eic_snr < settings$minimum_eic_snr) {
    result$status <- sprintf(
      "EIC SNR %.2f is below %.2f",
      result$eic_snr,
      settings$minimum_eic_snr
    )
    return(list(result = result, diagnostic = diagnostic))
  }

  mass_data <- collect_mass_points(
    ms_connection,
    metadata,
    ms1_rows,
    apex_rt_seconds,
    mz_range,
    settings$neighboring_scans_each_side,
    settings$minimum_eic_intensity
  )
  diagnostic$mass_points <- mass_data$points
  diagnostic$apex_spectrum <- mass_data$apex_spectrum

  required_mass_points <- 1L + 2L * settings$neighboring_scans_each_side
  if (nrow(mass_data$points) < required_mass_points) {
    result$apex_retention_time_min <- NA_real_
    result$status <- paste0(
      "fewer than ",
      required_mass_points,
      " contiguous scans contained a qualifying mass peak"
    )
    return(list(result = result, diagnostic = diagnostic))
  }

  result$observed_mz <- stats::weighted.mean(
    mass_data$points$mz,
    mass_data$points$intensity
  )
  result$absolute_ppm_difference <- absolute_ppm_error(
    theoretical_mz,
    result$observed_mz
  )
  result$accepted <- TRUE
  result$status <- "accepted"

  for (used_rt in mass_data$points$retention_time_min) {
    eic_row <- which.min(
      abs(diagnostic$eic$retention_time_min - used_rt)
    )
    diagnostic$eic$used_for_mass_average[eic_row] <- TRUE
  }

  list(result = result, diagnostic = diagnostic)
}

# Diagnostic plots

format_diagnostic_number <- function(value, digits = 2L) {
  if (length(value) == 0L || is.na(value)) {
    return("NA")
  }
  format(round(value, digits), nsmall = digits)
}

draw_empty_panel <- function(title, message) {
  graphics::plot.new()
  graphics::title(main = title)
  graphics::text(0.5, 0.5, message)
}

has_qualifying_apex_spectrum <- function(diagnostic) {
  spectrum <- diagnostic$apex_spectrum
  !is.null(spectrum) && nrow(spectrum) > 0L
}

write_diagnostic_report <- function(file, results, diagnostics, output_file,
                                    settings) {
  grDevices::pdf(output_file, width = 10, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (index in seq_len(nrow(results))) {
    result <- results[index, , drop = FALSE]
    diagnostic <- diagnostics[[index]]
    graphics::par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2, 0))

    eic <- diagnostic$eic
    if (nrow(eic) > 0L) {
      graphics::plot(
        eic$retention_time_min,
        eic$intensity,
        type = "l",
        xlab = "Retention time (min)",
        ylab = "EIC intensity",
        ylim = c(0, max(1, eic$intensity, na.rm = TRUE)),
        main = sprintf(
          "EIC %.4f m/z | SNR %s | %s",
          result$theoretical_mz,
          format_diagnostic_number(result$eic_snr),
          result$status
        )
      )
      if (is.finite(result$expected_rt_min)) {
        graphics::abline(
          v = result$expected_rt_min,
          col = "steelblue4",
          lty = 2
        )
      }
      if (is.finite(result$apex_retention_time_min)) {
        graphics::abline(
          v = result$apex_retention_time_min,
          col = "firebrick3",
          lty = 2
        )
      }
      used <- which(eic$used_for_mass_average)
      if (length(used) > 0L) {
        graphics::points(
          eic$retention_time_min[used],
          eic$intensity[used],
          pch = 19,
          col = "darkorange3"
        )
      }
    } else {
      draw_empty_panel(
        sprintf("EIC %.4f m/z", result$theoretical_mz),
        result$status
      )
    }

    spectrum <- diagnostic$apex_spectrum
    if (has_qualifying_apex_spectrum(diagnostic)) {
      plot_tolerance <- result$theoretical_mz *
        settings$spectrum_plot_ppm * 1e-6
      plotted <- spectrum[
        spectrum[, "mz"] >= result$theoretical_mz - plot_tolerance &
          spectrum[, "mz"] <= result$theoretical_mz + plot_tolerance, ,
        drop = FALSE
      ]
      if (nrow(plotted) > 0L) {
        graphics::plot(
          plotted[, "mz"],
          plotted[, "intensity"],
          type = "h",
          xlab = "m/z",
          ylab = "Intensity",
          ylim = c(0, max(1, plotted[, "intensity"], na.rm = TRUE)),
          main = paste0(
            "Apex MS1 spectrum | observed m/z ",
            format_diagnostic_number(result$observed_mz, 6L),
            " | absolute error ",
            format_diagnostic_number(
              result$absolute_ppm_difference,
              3L
            ),
            " ppm"
          )
        )
        graphics::abline(
          v = result$theoretical_mz,
          col = "steelblue4",
          lty = 2
        )
        if (is.finite(result$observed_mz)) {
          graphics::abline(
            v = result$observed_mz,
            col = "darkorange3",
            lwd = 2,
            lty = 2
          )
        }
      } else {
        draw_empty_panel(
          "Apex MS1 spectrum",
          "No signals in the plotting window"
        )
      }
    } else {
      draw_empty_panel(
        "Apex MS1 spectrum",
        "No qualifying apex spectrum"
      )
    }

    graphics::mtext(
      paste(basename(file), result$polarity, sep = " | "),
      side = 3,
      outer = TRUE,
      line = 0.5
    )
  }

  invisible(output_file)
}

# File calibration

calibrate_mzml_file <- function(file, output_directory,
                                positive_standards, negative_standards,
                                settings) {
  message("Calibrating: ", basename(file))
  ms_connection <- mzR::openMSfile(file)
  on.exit(mzR::close(ms_connection), add = TRUE)
  metadata <- mzR::header(ms_connection)

  ms1_rows <- which(
    metadata$msLevel == 1L &
      metadata$peaksCount > 0L &
      is.finite(metadata$retentionTime)
  )
  if (length(ms1_rows) == 0L) {
    stop("No non-empty MS1 spectra were found in ", basename(file), ".")
  }

  polarity <- detect_polarity(metadata, file)
  standards <- if (polarity == "positive") {
    positive_standards
  } else {
    negative_standards
  }
  full_rt_range <- range(metadata$retentionTime[ms1_rows])

  ms_data <- MSnbase::readMSData(
    files = file,
    mode = "onDisk",
    msLevel. = 1L,
    verbose = FALSE
  )

  measurements <- vector("list", nrow(standards))
  for (standard_index in seq_len(nrow(standards))) {
    measurements[[standard_index]] <- measure_standard(
      file = file,
      standard = standards[standard_index, , drop = FALSE],
      polarity = polarity,
      ms_data = ms_data,
      ms_connection = ms_connection,
      metadata = metadata,
      ms1_rows = ms1_rows,
      full_rt_range = full_rt_range,
      settings = settings
    )
  }

  results <- do.call(rbind, lapply(measurements, `[[`, "result"))
  diagnostics <- lapply(measurements, `[[`, "diagnostic")
  has_apex_spectrum <- vapply(
    diagnostics,
    has_qualifying_apex_spectrum,
    logical(1)
  )
  results$apex_retention_time_min[!has_apex_spectrum] <- NA_real_
  accepted <- results$accepted & is.finite(results$absolute_ppm_difference)
  file_average <- if (any(accepted)) {
    mean(results$absolute_ppm_difference[accepted])
  } else {
    NA_real_
  }
  results$file_average_absolute_ppm_difference <- file_average

  result_directory <- result_directory_for_file(file, output_directory)
  file_stem <- tools::file_path_sans_ext(basename(file))
  csv_file <- file.path(
    result_directory,
    paste0(file_stem, "_mass_accuracy.csv")
  )
  pdf_file <- file.path(
    result_directory,
    paste0(file_stem, "_mass_calibration_diagnostics.pdf")
  )

  csv_columns <- c(
    "theoretical_mz",
    "expected_rt_min",
    "observed_mz",
    "absolute_ppm_difference",
    "apex_retention_time_min",
    "file_average_absolute_ppm_difference"
  )
  utils::write.csv(
    results[, csv_columns, drop = FALSE],
    csv_file,
    row.names = FALSE,
    na = ""
  )

  if (settings$generate_diagnostic_pdf) {
    write_diagnostic_report(
      file,
      results,
      diagnostics,
      pdf_file,
      settings
    )
  }

  message(
    "  Polarity: ", polarity,
    " | accepted standards: ", sum(accepted), "/", nrow(results),
    " | mean absolute error: ",
    if (is.finite(file_average)) {
      paste0(format(round(file_average, 4), nsmall = 4), " ppm")
    } else {
      "not available"
    }
  )
  message("  CSV: ", csv_file)
  if (settings$generate_diagnostic_pdf) message("  PDF: ", pdf_file)

  diagnostic_pdf <- if (settings$generate_diagnostic_pdf) {
    pdf_file
  } else {
    NA_character_
  }
  data.frame(
    file = basename(file),
    polarity = polarity,
    accepted_standards = sum(accepted),
    total_standards = nrow(results),
    average_absolute_ppm_difference = file_average,
    csv_file = csv_file,
    diagnostic_pdf = diagnostic_pdf,
    stringsAsFactors = FALSE
  )
}

run_mass_calibration <- function(input_path, output_directory = NULL,
                                 positive_standards, negative_standards,
                                 settings) {
  check_required_packages()
  validate_configuration(
    input_path,
    output_directory,
    positive_standards,
    negative_standards,
    settings
  )
  files <- find_mzml_files(input_path)
  message("Found ", length(files), " mzML file(s). Processing sequentially.")

  summaries <- vector("list", length(files))
  for (file_index in seq_along(files)) {
    summaries[[file_index]] <- calibrate_mzml_file(
      files[file_index],
      output_directory,
      positive_standards,
      negative_standards,
      settings
    )
  }

  summary <- do.call(rbind, summaries)
  message("Mass calibration complete.")
  invisible(summary)
}

# Run

calibration_settings <- list(
  eic_search_ppm = eic_search_ppm,
  retention_time_tolerance_min = retention_time_tolerance_min,
  minimum_eic_intensity = minimum_eic_intensity,
  minimum_eic_scans = minimum_eic_scans,
  minimum_eic_snr = minimum_eic_snr,
  expected_peak_width_seconds = expected_peak_width_seconds,
  neighboring_scans_each_side = neighboring_scans_each_side,
  generate_diagnostic_pdf = generate_diagnostic_pdf,
  spectrum_plot_ppm = spectrum_plot_ppm
)

if (identical(environment(), globalenv())) {
  calibration_summary <- run_mass_calibration(
    input_path = input_path,
    output_directory = output_directory,
    positive_standards = positive_standards,
    negative_standards = negative_standards,
    settings = calibration_settings
  )
}
