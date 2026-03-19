#' Create a randomised wet-lab batch of samples and generate a PDF pick list
#'
#' Selects a balanced random batch of eligible samples for a chosen omic assay
#' (e.g. Proteomic), optionally matches to a biobank inventory to generate a tube
#' pick list, updates processing status in the sample table, and renders a PDF report.
#'
#' This function is designed to aid with wet-lab batching:
#' \itemize{
#'   \item Eligibility is controlled by a per-omic status column (e.g. Proteomic == "Yes")
#'   \item Randomisation can be balanced across one or more variables of interest
#'   \item Biobank lookup can skip used tubes and pick the next available replicate(s)
#' }
#'
#' @param samples A `data.frame`/`tibble` or a file path to `.csv`, `.tsv`, or `.xlsx`.
#'
#' @param id_col Name of the unique sample identifier column in `samples`.
#'   Defaults to `"Unique_ID"`. Values must be unique and non-missing.
#'
#' @param strata_cols Character vector of columns used to balance randomisation.
#'   Balancing is performed on weightings defined by these columns.
#'
#' @param omic Name of the omic status column in `samples` (e.g. "Proteomic").
#'   Used to filter eligible samples (default "Yes") and optionally update processed samples to no longer be eligible in future batches (default "No"). These strings can be edited via eligible_sample or not_eligible_sample arguments.
#'
#' @param batch_size Number of samples to select.
#'
#' @param seed Optional integer seed for reproducibility. If `NULL`, a seed is generated and returned.
#'
#' @param eligible_sample Character string to define eligibility for inclusion in batch. Default "Yes".
#'
#' @param not_eligible_sample Character string to define non-eligibility for inclusion in batch. Default "No".
#'
#' @param balance_strategy Either `"proportional"` (default) or `"stratified"`.
#'   `"stratified"` Allocates per-level targets for each balancing variable based on eligible cohort proportions, then selects a random set of samples that best matches these targets across all variables simultaneously (prioritising under-represented groups when trade-offs are required).
#'   `"proportional"` batches match the eligible cohort composition on average.
#'
#' @param sample_display_cols Optional columns to print in the PDF for the selected batch.
#'   If `NULL`, the report prints `id_col` plus `strata_cols`.
#'
#' @param biobank Optional `data.frame` or file path. If `NULL`, tube lookup is skipped
#'   and the report omits the biobank table.
#'
#' @param biobank_id_col Column in `biobank` containing the same unique IDs as `samples[[id_col]]`.
#'   Defaults to `"Unique_ID"`.
#'
#' @param biobank_omic_col Column in `biobank` indicating the omic the sample relates to (default "Application").
#'
#' @param tube_id_col Column in `biobank` giving the unique tube ID (default "Tube_Specific_ID").
#'
#' @param used_col Optional column in `biobank` indicating if tube has been used (default `NULL`).
#'   When supplied, only rows where `used_col == used_no` are considered available and selected tubes are updated to `used_yes`.
#'
#' @param used_yes Value in `used_col` column indicating if tube has been used (default "Yes").
#'   Only used when `used_col` is supplied.
#'
#' @param used_no Value in `used_col` column indicating available tubes (default "No").
#'   Only used when `used_col` is supplied.
#'
#' @param replicates_col Optional column indicating replicate number (default `NULL`).
#'   When supplied, the lowest replicate(s) are chosen first.
#'
#' @param tubes_per_sample Integer number of tubes to pick per sample (default 1).
#'   When `biobank` is provided, samples with fewer than this number of available tubes are reported in `insufficient_tubes` and are not marked as processed.
#'
#' @param position_cols Character vector of biobank location columns to include in the tube pick list (i.e. row, column, box, etc).
#'
#' @param pdf_file Output path for the PDF. If `NULL`, a timestamped file is created in the working directory.
#'   Relative paths are resolved from the caller's working directory.
#'
#' @param write_updates Logical; if `TRUE`, writes updated samples/biobank tables to disk (default `FALSE`).
#' omic col in samples file will be updated, and Batch_sample_processed will be updated (or created if not present) for successfully picked tubes. If `used_col` is supplied, it will also be updated in the biobank file. Files are only written after the PDF report renders successfully.
#'
#' @param update_in_place Logical; if `TRUE`, overwrite the original input files (only when inputs are file paths).
#'   Default `FALSE`, which writes new `_UPDATED` files instead.
#'
#' @param samples_out Optional output path for updated samples table. If `NULL` and `write_updates=TRUE`,
#'   an `_UPDATED` path is generated (unless `update_in_place=TRUE`).
#'
#' @param biobank_out Optional output path for updated biobank table. Same rules as `samples_out`.
#'
#' @param quiet Logical; reduce rmarkdown output (default TRUE).
#'
#' @return A list with elements:
#' \describe{
#'   \item{batch}{Selected rows from the eligible subset of `samples`.}
#'   \item{updated_samples}{samples table with omic status updated for processed samples.}
#'   \item{biobank_picklist}{Tube pick list for samples with enough available tubes (or `NULL` if `biobank` not provided).}
#'   \item{updated_biobank}{Biobank table with selected tubes annotated for the batch and, when `used_col` is supplied, updated to `used_yes` (or `NULL`).}
#'   \item{unmatched_samples}{Selected samples that could not be matched to any available tube (or `NULL`).}
#'   \item{insufficient_tubes}{Selected samples with fewer than `tubes_per_sample` available tubes (or `NULL`).}
#'   \item{seed}{The seed used.}
#'   \item{pdf}{Path to the rendered PDF (or `NA` if PDF rendering dependencies are missing).}
#'   \item{written_files}{Named character vector of written files (empty if `write_updates=FALSE`).}
#' }
#'
#' @importFrom stats setNames
#' @importFrom utils head
#'
#' @examples
#' samples_path <- system.file("extdata", "samples.csv", package = "batchcreator")
#' biobank_path <- system.file("extdata", "biobank.csv", package = "batchcreator")
#'
#' \dontrun{
#' out <- create_batch(
#'   samples = samples_path,
#'   strata_cols = c("Sample_group", "Sex"),
#'   omic = "Proteomic",
#'   batch_size = 10,
#'   biobank = biobank_path,
#'   used_col = "Used",
#'   replicates_col = "Replicates",
#'   pdf_file = tempfile(fileext = ".pdf")
#' )
#' }
#'
#' @export
create_batch <- function(
    samples,
    id_col = "Unique_ID",
    strata_cols,
    omic,
    batch_size,
    seed = NULL,
    eligible_sample = "Yes",
    not_eligible_sample = "No",
    balance_strategy = c("proportional", "stratified"),
    sample_display_cols = NULL,
    biobank = NULL,
    biobank_id_col = "Unique_ID",
    biobank_omic_col = "Application",
    tube_id_col = "Tube_Specific_ID",
    used_col = NULL,
    used_yes = "Yes",
    used_no = "No",
    replicates_col = NULL,
    tubes_per_sample = 1,
    position_cols = c("Freezer", "Stage", "Rack", "Drawer", "Slot", "Box", "Position"),
    pdf_file = NULL,
    write_updates = FALSE,
    update_in_place = FALSE,
    samples_out = NULL,
    biobank_out = NULL,
    quiet = TRUE
) {
  balance_strategy <- match.arg(balance_strategy)

  if (missing(strata_cols) || length(strata_cols) < 1L) {
    stop("`strata_cols` must be provided to define how batches are randomised (e.g. c('Group','Sex','Replicate','Cohort', etc)).")
  }
  if (missing(omic) || length(omic) != 1L) {
    stop("`omic` must be a single column name (e.g. 'Proteomic').")
  }

  samples_is_path <- is.character(samples) && length(samples) == 1L
  biobank_is_path  <- is.character(biobank)  && length(biobank)  == 1L

  samples_df <- .read_tabular(samples, object_name = "samples")
  .assert_has_cols(samples_df, c(id_col, strata_cols, omic), object_name = "samples")

  .assert_unique_id(samples_df, id_col, object_name = "samples")

  # Validate omic status values
  if (!all(samples_df[[omic]] %in% c(eligible_sample, not_eligible_sample))) {
    bad <- unique(samples_df[[omic]][!samples_df[[omic]] %in% c(eligible_sample, not_eligible_sample)])
    stop(
      "Unexpected values in samples[[", omic, "]]. Expected only ",
      shQuote(eligible_sample), " / ", shQuote(not_eligible_sample),
      ". Found: ", paste(shQuote(bad), collapse = ", "),
      ". If your file uses different labels, set arguments `eligible_sample` / `not_eligible_sample`."
    )
  }

  eligible <- samples_df[samples_df[[omic]] == eligible_sample, , drop = FALSE]
  if (nrow(eligible) == 0L) {
    stop("No eligible samples found (", omic, " == ", shQuote(eligible_sample), "). Nothing to batch.")
  }

  # Batch size checks
  if (missing(batch_size)) {
    stop("`batch_size` must be provided.")
  }
  batch_size <- .validate_positive_whole_number(batch_size, "batch_size")

  note <- NULL
  if (nrow(eligible) < batch_size) {
    note <- paste0(
      "Requested batch_size=", batch_size,
      " exceeds eligible population (n=", nrow(eligible),
      "). Using batch_size=", nrow(eligible), " (all eligible samples)."
    )
    batch_size <- nrow(eligible)
  }

  # Seed
  if (is.null(seed)) {
    seed <- sample.int(100000L, 1L)
  } else {
    seed <- .validate_seed(seed)
  }
  set.seed(seed)

  # Decide PDF filename early (needed for biobank annotation)
  if (is.null(pdf_file)) {
    stamp <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
    pdf_file <- file.path(getwd(), paste0("Batch_", omic, "_", stamp, ".pdf"))
  } else if (!.is_absolute_path(pdf_file)) {
    pdf_file <- file.path(getwd(), pdf_file)
  }
  batch_name <- basename(pdf_file)

  # Sampling
  idx <- switch(
    balance_strategy,
    proportional = sample.int(nrow(eligible), size = batch_size, replace = FALSE),
    stratified   = .marginal_balanced_sample_idx(eligible, strata_cols, batch_size)
  )

  batch <- eligible[idx, , drop = FALSE]

  processed_ids <- batch[[id_col]]
  selected_tubes <- character(0)

  # Optional biobank lookup + Used updates
  biobank_df <- NULL
  biobank_picklist <- NULL
  updated_biobank <- NULL
  unmatched_samples <- NULL
  insufficient_tubes <- NULL

  if (!is.null(biobank)) {
    biobank_df <- .read_tabular(biobank, object_name = "biobank")
    .assert_has_cols(
      biobank_df,
      unique(c(biobank_id_col, biobank_omic_col, tube_id_col, position_cols)),
      object_name = "biobank"
    )
    if (!is.null(used_col)) {
      .assert_has_cols(biobank_df, used_col, object_name = "biobank")
    }
    if (!is.null(replicates_col)) {
      .assert_has_cols(biobank_df, replicates_col, object_name = "biobank")
    }

    # ID sanity (doesn't have to be unique in biobank because there are replicates)
    if (.has_missing_or_empty(biobank_df[[biobank_id_col]])) {
      stop("`biobank[[", biobank_id_col, "]]` contains missing/empty IDs. Please fix to ensure no empty/missing values")
    }
    if (.has_missing_or_empty(biobank_df[[tube_id_col]])) {
      stop("`biobank[[", tube_id_col, "]]` contains missing/empty tube IDs. Please fix to ensure no empty/missing values")
    }
    if (anyDuplicated(biobank_df[[tube_id_col]]) > 0L) {
      stop(
        "`biobank[[", tube_id_col, "]]` must be unique. ",
        "Tube IDs are used to mark selected inventory rows as used, so duplicates would update multiple rows."
      )
    }

    processed_ids <- batch[[id_col]][0]

    # Filter to omic/Omic
    biobank_sub <- biobank_df[biobank_df[[biobank_omic_col]] == omic, , drop = FALSE]
    if (nrow(biobank_sub) == 0L) {
      warning("Biobank provided, but no rows matched ", biobank_omic_col, " == ", shQuote(omic), ".")
      unmatched_samples <- batch[, id_col, drop = FALSE]
    } else {
      # Keep only available tubes if Used present
      sub_available <- biobank_sub
      if (!is.null(used_col) && used_col %in% names(sub_available)) {
        sub_available <- sub_available[sub_available[[used_col]] == used_no, , drop = FALSE]
      }

      # Merge selected samples with available tube candidates
      merged <- merge(
        batch,
        sub_available,
        by.x = id_col,
        by.y = biobank_id_col,
        all.x = TRUE,
        sort = FALSE
      )

      # Samples with zero available tubes
      no_tube <- is.na(merged[[tube_id_col]])
      if (any(no_tube)) {
        unmatched_samples <- unique(merged[no_tube, id_col, drop = FALSE])
      }
      matched <- merged[!no_tube, , drop = FALSE]

      tubes_per_sample <- .validate_positive_whole_number(tubes_per_sample, "tubes_per_sample")

      # Order candidates: by replicate if available, otherwise stable by tube id
      if (!is.null(replicates_col) && replicates_col %in% names(matched)) {
        matched <- matched[order(matched[[id_col]], matched[[replicates_col]]), , drop = FALSE]
      } else {
        matched <- matched[order(matched[[id_col]], matched[[tube_id_col]]), , drop = FALSE]
      }

      # Check insufficient tubes (sample has < tubes_per_sample available)
      counts <- table(matched[[id_col]])
      insufficient_ids <- names(counts)[counts < tubes_per_sample]
      sufficient_ids <- names(counts)[counts >= tubes_per_sample]
      if (length(insufficient_ids) > 0L) {
        insufficient_tubes <- data.frame(
          Unique_ID = insufficient_ids,
          available_tubes = as.integer(counts[insufficient_ids]),
          required_tubes = tubes_per_sample,
          stringsAsFactors = FALSE
        )
        names(insufficient_tubes)[1] <- id_col
      }

      if (length(sufficient_ids) > 0L) {
        matched_sufficient <- matched[matched[[id_col]] %in% sufficient_ids, , drop = FALSE]
        pick <- .take_first_n_per_id(matched_sufficient, id_col = id_col, n = tubes_per_sample)
        processed_ids <- unique(pick[[id_col]])
        selected_tubes <- unique(pick[[tube_id_col]])
      } else {
        pick <- matched[FALSE, , drop = FALSE]
      }

      # Tidy picklist columns
      keep_cols <- unique(c(id_col, strata_cols, tube_id_col, replicates_col, position_cols))
      keep_cols <- keep_cols[keep_cols %in% names(pick)]
      biobank_picklist <- pick[, keep_cols, drop = FALSE]
    }
  }

  # Render PDF report

  pdf_path <- .render_batch_report(
    pdf_file = pdf_file,
    omic = omic,
    seed = seed,
    note = note,
    eligible = eligible,
    batch = batch,
    id_col = id_col,
    strata_cols = strata_cols,
    sample_display_cols = sample_display_cols,
    biobank_picklist = biobank_picklist,
    unmatched_samples = unmatched_samples,
    insufficient_tubes = insufficient_tubes,
    quiet = quiet
  )

  report_rendered <- !is.na(pdf_path) && nzchar(pdf_path)

  # Update samples status
  updated_samples <- samples_df
  if (length(processed_ids) > 0L) {
    sel <- updated_samples[[id_col]] %in% processed_ids
    updated_samples[sel, omic] <- not_eligible_sample
  }

  # Prepare updated biobank after report status is known
  if (!is.null(biobank_df)) {
    updated_biobank <- biobank_df

    if (length(selected_tubes) > 0L) {
      to_update <- updated_biobank[[tube_id_col]] %in% selected_tubes

      if (!is.null(used_col) && used_col %in% names(updated_biobank)) {
        updated_biobank[to_update, used_col] <- used_yes
      }

      if (isTRUE(report_rendered)) {
        if (!"Batch_sample_processed" %in% names(updated_biobank)) {
          updated_biobank[["Batch_sample_processed"]] <- NA_character_
        }
        updated_biobank[to_update, "Batch_sample_processed"] <- batch_name
      }
    }
  }

  # Decide output paths + write updates if requested
  written_files <- character(0)

  if (isTRUE(write_updates) && !isTRUE(report_rendered)) {
    warning("Skipping file writes because the PDF report was not rendered successfully. No files were modified on disk.")
  }

  if (isTRUE(write_updates) && isTRUE(report_rendered)) {
    if (samples_is_path) {
      if (is.null(samples_out)) {
        samples_out <- if (isTRUE(update_in_place)) {
          samples
        } else {
          .with_suffix(samples, "_UPDATED")
        }
      }
      .write_tabular(updated_samples, samples_out, object_name = "updated samples")
      written_files["samples"] <- samples_out
    } else if (!is.null(samples_out)) {
      .write_tabular(updated_samples, samples_out, object_name = "updated samples")
      written_files["samples"] <- samples_out
    }

    if (!is.null(biobank) && !is.null(updated_biobank)) {
      if (biobank_is_path) {
        if (is.null(biobank_out)) {
          biobank_out <- if (isTRUE(update_in_place)) {
            biobank
          } else {
            .with_suffix(biobank, "_UPDATED")
          }
        }
        .write_tabular(updated_biobank, biobank_out, object_name = "updated biobank")
        written_files["biobank"] <- biobank_out
      } else if (!is.null(biobank_out)) {
        .write_tabular(updated_biobank, biobank_out, object_name = "updated biobank")
        written_files["biobank"] <- biobank_out
      }
    }
  }

  invisible(list(
    batch = batch,
    updated_samples = updated_samples,
    biobank_picklist = biobank_picklist,
    updated_biobank = updated_biobank,
    unmatched_samples = unmatched_samples,
    insufficient_tubes = insufficient_tubes,
    seed = seed,
    pdf = pdf_path,
    written_files = written_files
  ))
}

# ---- internal helpers ----

.with_suffix <- function(path, suffix) {
  ext <- tools::file_ext(path)
  base <- sub(paste0("\\.", ext, "$"), "", path)
  paste0(base, suffix, ".", ext)
}

.is_absolute_path <- function(path) {
  grepl("^(?:[A-Za-z]:[\\\\/]|/|\\\\\\\\)", path)
}

.has_missing_or_empty <- function(x) {
  if (anyNA(x)) return(TRUE)
  any(trimws(as.character(x)) == "")
}

.validate_positive_whole_number <- function(x, arg_name) {
  tol <- sqrt(.Machine$double.eps)

  if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x)) {
    stop("`", arg_name, "` must be a single positive whole number.")
  }
  if (x <= 0 || abs(x - round(x)) > tol) {
    stop("`", arg_name, "` must be a single positive whole number.")
  }

  as.integer(round(x))
}

.validate_seed <- function(seed) {
  tol <- sqrt(.Machine$double.eps)

  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed) || !is.finite(seed)) {
    stop("`seed` must be NULL or a single whole number.")
  }
  if (abs(seed - round(seed)) > tol) {
    stop("`seed` must be NULL or a single whole number.")
  }

  as.integer(round(seed))
}

.assert_has_cols <- function(df, cols, object_name = "data") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0L) {
    stop(
      "`", object_name, "` is missing required column(s): ",
      paste(shQuote(missing), collapse = ", "),
      ".\nAvailable columns: ",
      paste(shQuote(names(df)), collapse = ", ")
    )
  }
}

.assert_unique_id <- function(df, id_col, object_name = "data") {
  if (.has_missing_or_empty(df[[id_col]])) {
    stop("`", object_name, "[[", id_col, "]]` contains missing/empty IDs. Please fix to ensure no empty/missing values")
  }
  if (anyDuplicated(df[[id_col]]) > 0L) {
    stop(
      "`", object_name, "[[", id_col, "]]` is not unique. ",
      "Create a single unique sample ID column (e.g. Unique_ID) and use that."
    )
  }
}

.read_tabular <- function(x, object_name = "data") {
  if (inherits(x, "data.frame")) return(x)

  if (!is.character(x) || length(x) != 1L) {
    stop("`", object_name, "` must be a data.frame or a single file path.")
  }
  if (!file.exists(x)) stop("File not found for `", object_name, "`: ", x)

  ext <- tolower(tools::file_ext(x))
  if (ext == "csv") {
    return(utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext %in% c("tsv", "txt")) {
    return(utils::read.delim(x, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Reading Excel requires the 'readxl' package. Install it or save the file as CSV/TSV.")
    }
    return(as.data.frame(readxl::read_excel(x)))
  }

  stop("Unsupported file extension for `", object_name, "`: .", ext,
       " (supported: csv, tsv/txt, xlsx/xls).")
}

.write_tabular <- function(df, path, object_name = "data") {
  ext <- tolower(tools::file_ext(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (ext == "csv") {
    utils::write.csv(df, path, row.names = FALSE)
    return(invisible(path))
  }
  if (ext %in% c("tsv", "txt")) {
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    return(invisible(path))
  }
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Writing Excel requires the 'openxlsx' package. Install it or write to CSV/TSV instead.")
    }
    openxlsx::write.xlsx(df, path, overwrite = TRUE)
    return(invisible(path))
  }

  stop("Unsupported output extension for `", object_name, "`: .", ext,
       " (supported: csv, tsv/txt, xlsx/xls).")
}

.make_strata_key <- function(df, strata_cols) {
  parts <- lapply(strata_cols, function(cc) as.character(df[[cc]]))
  do.call(paste, c(parts, sep = " | "))
}

.take_first_n_per_id <- function(df, id_col, n = 1L) {
  n <- as.integer(n)
  split_idx <- split(seq_len(nrow(df)), df[[id_col]])
  keep <- unlist(lapply(split_idx, function(ii) head(ii, n)), use.names = FALSE)
  df[keep, , drop = FALSE]
}

.render_batch_report <- function(
    pdf_file,
    omic,
    seed,
    note,
    eligible,
    batch,
    id_col,
    strata_cols,
    sample_display_cols,
    biobank_picklist,
    unmatched_samples,
    insufficient_tubes,
    quiet = TRUE
) {
  if (!requireNamespace("rmarkdown", quietly = TRUE) ||
      !requireNamespace("knitr", quietly = TRUE)) {
    warning("PDF report not rendered because 'rmarkdown' and/or 'knitr' are not installed. Returning data only.")
    return(NA_character_)
  }

  template <- system.file("templates", "batch_report.Rmd", package = "batchcreator")
  if (template == "") {
    local_template <- file.path("inst", "templates", "batch_report.Rmd")
    if (file.exists(local_template)) {
      template <- normalizePath(local_template, winslash = "/", mustWork = TRUE)
    }
  }
  if (template == "") {
    warning("Report template not found in package (inst/templates/batch_report.Rmd). Returning data only.")
    return(NA_character_)
  }

  eligible_summary <- .count_by_strata(eligible, strata_cols)
  batch_summary <- .count_by_strata(batch, strata_cols)

  if (is.null(sample_display_cols)) {
    sample_display_cols <- unique(c(id_col, strata_cols))
  }
  sample_display_cols <- sample_display_cols[sample_display_cols %in% names(batch)]
  batch_table <- batch[, sample_display_cols, drop = FALSE]

  e <- new.env(parent = baseenv())
  e$omic <- omic
  e$seed <- seed
  e$note <- note
  e$strata_cols <- strata_cols
  e$eligible <- eligible
  e$batch <- batch
  e$eligible_summary <- eligible_summary
  e$batch_summary <- batch_summary
  e$batch_table <- batch_table
  e$biobank_picklist <- biobank_picklist
  e$unmatched_samples <- unmatched_samples
  e$insufficient_tubes <- insufficient_tubes

  pdf_dir <- dirname(pdf_file)
  if (!nzchar(pdf_dir) || identical(pdf_dir, ".")) {
    pdf_dir <- getwd()
  }
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_dir <- normalizePath(pdf_dir, winslash = "/", mustWork = TRUE)
  pdf_name <- basename(pdf_file)

  rendered <- tryCatch(
    rmarkdown::render(
      input = template,
      output_file = pdf_name,
      output_dir = pdf_dir,
      quiet = quiet,
      envir = e
    ),
    error = function(e) {
      warning("PDF report was not rendered successfully: ", conditionMessage(e), ". Returning data only.")
      NA_character_
    }
  )

  if (is.na(rendered)) {
    return(NA_character_)
  }

  normalizePath(rendered, winslash = "/", mustWork = FALSE)
}

.count_by_strata <- function(df, strata_cols) {
  key <- .make_strata_key(df, strata_cols)
  out <- as.data.frame(table(key), stringsAsFactors = FALSE)
  names(out) <- c("strata", "n")

  parts <- strsplit(out$strata, " \\| ")
  parts <- do.call(rbind, lapply(parts, function(x) {
    length(x) <- length(strata_cols)
    x
  }))
  parts <- as.data.frame(parts, stringsAsFactors = FALSE)
  names(parts) <- strata_cols
  cbind(parts, n = out$n)
}

.marginal_balanced_sample_idx <- function(df, strata_cols, batch_size) {
  # Build per-sample membership keys (one per column), e.g. "Sex=F"
  key_mat <- sapply(strata_cols, function(cc) paste0(cc, "=", as.character(df[[cc]])))
  if (is.null(dim(key_mat))) key_mat <- matrix(key_mat, ncol = 1)

  n <- nrow(df)
  if (batch_size > n) stop("`batch_size` exceeds eligible population for stratified sampling.")

  # Targets per column-level (each column sums to batch_size)
  targets <- .marginal_targets(key_mat, batch_size)

  # Availability + weights (rare levels weighted higher)
  avail <- table(as.vector(key_mat))
  w <- 1 / pmax(as.numeric(avail), 1)
  names(w) <- names(avail)

  # Deficits: how many more of each level we still want
  deficit <- setNames(integer(length(avail)), names(avail))
  deficit[names(targets)] <- targets

  remaining <- seq_len(n)
  chosen <- integer(0)

  for (step in seq_len(batch_size)) {
    km <- key_mat[remaining, , drop = FALSE]

    # Marginal improvement from picking a sample = sum weights for unmet levels
    keys_vec <- as.vector(km)

    ww <- matrix(w[keys_vec], nrow = nrow(km), ncol = ncol(km))
    dd <- matrix(deficit[keys_vec] > 0, nrow = nrow(km), ncol = ncol(km))

    ww[is.na(ww)] <- 0
    dd[is.na(dd)] <- FALSE

    gain <- rowSums(ww * dd)

    if (max(gain) <= 0) {
      pick <- sample(remaining, 1) # no unmet deficits left (or not achievable)
    } else {
      cand <- remaining[gain == max(gain)]
      pick <- sample(cand, 1)
    }

    chosen <- c(chosen, pick)

    # Update deficits for this sample's levels (one per column)
    keys <- key_mat[pick, ]
    deficit[keys] <- deficit[keys] - 1L

    remaining <- remaining[remaining != pick]
  }

  chosen
}

.marginal_targets <- function(key_mat, batch_size) {
  # key_mat: n x K character matrix ("col=value")
  K <- ncol(key_mat)
  out <- integer(0)

  for (j in seq_len(K)) {
    counts <- table(key_mat[, j])
    props <- counts / sum(counts)
    expected <- as.numeric(props) * batch_size

    alloc <- .allocate_targets_lrm(
      expected = expected,
      avail = as.integer(counts),
      names = names(counts),
      total = batch_size
    )

    out <- c(out, alloc)
  }

  out
}

.allocate_targets_lrm <- function(expected, avail, names, total) {
  # Largest remainder with capacity + small-group tie-break
  base <- floor(expected)
  base <- pmin(base, avail)

  remaining <- total - sum(base)
  frac <- expected - floor(expected)

  # Distribute remaining to the largest fractional parts; if tied, prefer smaller groups
  ord <- order(-frac, avail, expected)  # small-group priority via avail/expected

  i <- 1L
  while (remaining > 0) {
    idx <- ord[i]
    if (base[idx] < avail[idx]) {
      base[idx] <- base[idx] + 1L
      remaining <- remaining - 1L
    }
    i <- i + 1L
    if (i > length(ord)) i <- 1L
    if (all(base >= avail)) break
  }

  stats::setNames(as.integer(base), names)
}

