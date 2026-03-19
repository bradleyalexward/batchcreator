#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 140)
set.seed(1)

OUTPUT_DIR <- file.path("inst", "scripts", "test_outputs", "create_batch")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (exists("create_batch", mode = "function")) {
  message("create_batch() already exists in session; reloading from R/create_batch.R")
}
source(file.path("R", "create_batch.R"))

samples0 <- read.csv("inst/extdata/samples.csv", check.names = FALSE, stringsAsFactors = FALSE)
biobank0 <- read.csv("inst/extdata/biobank.csv", check.names = FALSE, stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

check <- function(cond, fmt, ...) {
  if (!isTRUE(cond)) {
    stop(sprintf(fmt, ...), call. = FALSE)
  }
}

clean_text <- function(x, n = 1L) {
  if (length(x) == 0L) {
    return("")
  }
  paste(x[seq_len(min(length(x), n))], collapse = " | ")
}

same_df <- function(x, y) {
  identical(
    as.data.frame(x, stringsAsFactors = FALSE),
    as.data.frame(y, stringsAsFactors = FALSE)
  )
}

copy_inputs <- function(case_dir, samples_df = samples0, biobank_df = biobank0) {
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  samples_path <- file.path(case_dir, "samples.csv")
  write.csv(samples_df, samples_path, row.names = FALSE)

  biobank_path <- NULL
  if (!is.null(biobank_df)) {
    biobank_path <- file.path(case_dir, "biobank.csv")
    write.csv(biobank_df, biobank_path, row.names = FALSE)
  }

  list(
    samples = samples_path,
    biobank = biobank_path
  )
}

read_csv_safe <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

make_partial_inventory_data <- function() {
  ids <- samples0$Unique_ID[1:4]
  custom_samples <- samples0[match(ids, samples0$Unique_ID), , drop = FALSE]
  custom_samples$Proteomic <- "Yes"

  custom_biobank <- data.frame(
    Unique_ID = c(ids[1], ids[2], ids[2], ids[4], ids[4], ids[1]),
    Application = c("Proteomic", "Proteomic", "Proteomic", "Proteomic", "Proteomic", "Metabolomic"),
    Replicates = c(1, 1, 2, 1, 2, 1),
    Tube_Specific_ID = c("T1", "T2", "T3", "T4", "T5", "T6"),
    Freezer = c("F1", "F1", "F1", "F1", "F1", "F1"),
    Stage = c(1, 1, 1, 1, 1, 1),
    Rack = c(1, 1, 1, 1, 1, 1),
    Drawer = c(1, 1, 1, 1, 1, 1),
    Slot = c(1, 1, 1, 1, 1, 1),
    Box = c("Box1", "Box1", "Box1", "Box1", "Box1", "Box2"),
    Position = c("1A", "2A", "3A", "4A", "5A", "1B"),
    Used = c("No", "No", "No", "No", "No", "No"),
    Batch_sample_processed = rep("", 6),
    stringsAsFactors = FALSE
  )

  list(
    samples = custom_samples,
    biobank = custom_biobank,
    ids = ids
  )
}

cases <- list()
case_idx <- 0L

add_case <- function(name, runner, expect_error = FALSE, notes = "") {
  case_idx <<- case_idx + 1L
  cases[[case_idx]] <<- list(
    id = case_idx,
    name = name,
    runner = runner,
    expect_error = expect_error,
    notes = notes
  )
}

add_case(
  name = "baseline_proportional_biobank_dryrun",
  notes = "Baseline proportional run against extdata with biobank matching",
  runner = function(case_dir) {
    paths <- copy_inputs(case_dir)
    pdf_file <- file.path(case_dir, "baseline_proportional.pdf")

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 10,
      seed = 117,
      biobank = paths$biobank,
      biobank_id_col = "Unique_ID",
      biobank_omic_col = "Application",
      tube_id_col = "Tube_Specific_ID",
      used_col = "Used",
      used_yes = "Yes",
      used_no = "No",
      replicates_col = "Replicates",
      tubes_per_sample = 1,
      position_cols = c("Freezer", "Stage", "Rack", "Drawer", "Slot", "Box", "Position"),
      pdf_file = pdf_file,
      write_updates = FALSE,
      quiet = TRUE
    )

    check(is.list(res), "Expected list result from create_batch().")
    check(nrow(res$batch) == 10L, "Expected batch size 10, got %d.", nrow(res$batch))
    check(!is.null(res$biobank_picklist), "Expected a non-NULL biobank picklist.")
    check(all(res$biobank_picklist$Unique_ID %in% res$batch$Unique_ID), "Picklist contains IDs outside the selected batch.")

    processed_ids <- unique(res$biobank_picklist$Unique_ID)
    if (length(processed_ids) > 0L) {
      check(
        all(res$updated_samples[res$updated_samples$Unique_ID %in% processed_ids, "Proteomic"] == "No"),
        "Processed IDs were not marked as No in updated_samples."
      )
    }

    held_ids <- unique(c(
      res$unmatched_samples$Unique_ID %||% character(0),
      res$insufficient_tubes$Unique_ID %||% character(0)
    ))
    if (length(held_ids) > 0L) {
      check(
        all(res$updated_samples[res$updated_samples$Unique_ID %in% held_ids, "Proteomic"] == "Yes"),
        "Unmatched or insufficient IDs should remain eligible in updated_samples."
      )
    }

    list(
      n_batch = nrow(res$batch),
      n_picklist = nrow(res$biobank_picklist),
      n_unmatched = nrow(res$unmatched_samples %||% data.frame()),
      n_insufficient = nrow(res$insufficient_tubes %||% data.frame()),
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = FALSE,
      wrote_biobank = FALSE
    )
  }
)

add_case(
  name = "baseline_stratified_biobank_dryrun",
  notes = "Baseline stratified run against extdata",
  runner = function(case_dir) {
    paths <- copy_inputs(case_dir)
    pdf_file <- file.path(case_dir, "baseline_stratified.pdf")

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 10,
      seed = 118,
      balance_strategy = "stratified",
      biobank = paths$biobank,
      pdf_file = pdf_file,
      write_updates = FALSE,
      quiet = TRUE
    )

    check(nrow(res$batch) == 10L, "Expected stratified batch size 10, got %d.", nrow(res$batch))
    check(!is.null(res$biobank_picklist), "Expected a non-NULL biobank picklist for stratified case.")

    list(
      n_batch = nrow(res$batch),
      n_picklist = nrow(res$biobank_picklist),
      n_unmatched = nrow(res$unmatched_samples %||% data.frame()),
      n_insufficient = nrow(res$insufficient_tubes %||% data.frame()),
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = FALSE,
      wrote_biobank = FALSE
    )
  }
)

add_case(
  name = "no_biobank_dryrun",
  notes = "Selected samples should still be marked in updated_samples when biobank lookup is skipped",
  runner = function(case_dir) {
    paths <- copy_inputs(case_dir, biobank_df = NULL)
    pdf_file <- file.path(case_dir, "no_biobank.pdf")

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 5,
      seed = 119,
      biobank = NULL,
      pdf_file = pdf_file,
      write_updates = FALSE,
      quiet = TRUE
    )

    selected_ids <- res$batch$Unique_ID
    check(is.null(res$biobank_picklist), "Expected NULL biobank picklist when biobank is omitted.")
    check(length(res$written_files) == 0L, "Expected no written files in dry-run mode.")
    check(
      all(res$updated_samples[res$updated_samples$Unique_ID %in% selected_ids, "Proteomic"] == "No"),
      "Selected samples should be marked as No when no biobank is supplied."
    )

    list(
      n_batch = nrow(res$batch),
      n_picklist = 0L,
      n_unmatched = 0L,
      n_insufficient = 0L,
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = FALSE,
      wrote_biobank = FALSE
    )
  }
)

add_case(
  name = "write_updates_only_after_successful_pdf",
  notes = "When PDF rendering fails, no files should be modified on disk",
  runner = function(case_dir) {
    paths <- copy_inputs(case_dir)
    samples_before <- read_csv_safe(paths$samples)
    biobank_before <- read_csv_safe(paths$biobank)
    pdf_file <- file.path(case_dir, "write_updates.pdf")
    samples_out <- file.path(case_dir, "samples_UPDATED.csv")
    biobank_out <- file.path(case_dir, "biobank_UPDATED.csv")

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 5,
      seed = 120,
      biobank = paths$biobank,
      used_col = "Used",
      replicates_col = "Replicates",
      pdf_file = pdf_file,
      write_updates = TRUE,
      update_in_place = FALSE,
      samples_out = samples_out,
      biobank_out = biobank_out,
      quiet = TRUE
    )

    if (!is.na(res$pdf)) {
      check(file.exists(samples_out), "Expected updated samples file to be written when PDF exists.")
      check(file.exists(biobank_out), "Expected updated biobank file to be written when PDF exists.")

      samples_written <- read_csv_safe(samples_out)
      biobank_written <- read_csv_safe(biobank_out)

      check(same_df(samples_written, res$updated_samples), "Written samples file does not match returned updated_samples.")
      check(same_df(biobank_written, res$updated_biobank), "Written biobank file does not match returned updated_biobank.")
      check(length(res$written_files) == 2L, "Expected two written files when PDF exists.")
    } else {
      samples_after <- read_csv_safe(paths$samples)
      biobank_after <- read_csv_safe(paths$biobank)

      check(length(res$written_files) == 0L, "Expected no written_files entries when PDF render fails.")
      check(!file.exists(samples_out), "Did not expect samples_out to exist when PDF render fails.")
      check(!file.exists(biobank_out), "Did not expect biobank_out to exist when PDF render fails.")
      check(same_df(samples_before, samples_after), "Samples input should remain unchanged when PDF render fails.")
      check(same_df(biobank_before, biobank_after), "Biobank input should remain unchanged when PDF render fails.")
    }

    list(
      n_batch = nrow(res$batch),
      n_picklist = nrow(res$biobank_picklist %||% data.frame()),
      n_unmatched = nrow(res$unmatched_samples %||% data.frame()),
      n_insufficient = nrow(res$insufficient_tubes %||% data.frame()),
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = file.exists(samples_out),
      wrote_biobank = file.exists(biobank_out)
    )
  }
)

add_case(
  name = "oversize_batch_reduces_to_all_eligible",
  notes = "Requested batch_size above eligible count should return all eligible samples",
  runner = function(case_dir) {
    paths <- copy_inputs(case_dir, biobank_df = NULL)
    pdf_file <- file.path(case_dir, "oversize_batch.pdf")
    eligible_n <- sum(samples0$Proteomic == "Yes", na.rm = TRUE)

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = eligible_n + 100L,
      seed = 121,
      biobank = NULL,
      pdf_file = pdf_file,
      write_updates = FALSE,
      quiet = TRUE
    )

    check(nrow(res$batch) == eligible_n, "Expected %d eligible samples in oversize batch, got %d.", eligible_n, nrow(res$batch))

    list(
      n_batch = nrow(res$batch),
      n_picklist = 0L,
      n_unmatched = 0L,
      n_insufficient = 0L,
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = FALSE,
      wrote_biobank = FALSE
    )
  }
)

add_case(
  name = "fractional_batch_size_errors",
  expect_error = TRUE,
  notes = "batch_size should reject fractional values",
  runner = function(case_dir) {
    paths <- copy_inputs(case_dir, biobank_df = NULL)
    create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 0.5,
      biobank = NULL,
      quiet = TRUE
    )
  }
)

add_case(
  name = "fractional_tubes_per_sample_errors",
  expect_error = TRUE,
  notes = "tubes_per_sample should reject fractional values",
  runner = function(case_dir) {
    paths <- copy_inputs(case_dir)
    create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 2,
      seed = 122,
      biobank = paths$biobank,
      tubes_per_sample = 1.5,
      quiet = TRUE
    )
  }
)

add_case(
  name = "duplicate_sample_id_errors",
  expect_error = TRUE,
  notes = "Duplicate sample IDs should be rejected",
  runner = function(case_dir) {
    dup_samples <- samples0
    dup_samples$Unique_ID[2] <- dup_samples$Unique_ID[1]
    paths <- copy_inputs(case_dir, samples_df = dup_samples, biobank_df = NULL)
    create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 1,
      biobank = NULL,
      quiet = TRUE
    )
  }
)

add_case(
  name = "duplicate_tube_id_errors",
  expect_error = TRUE,
  notes = "Duplicate biobank tube IDs should be rejected before updates are prepared",
  runner = function(case_dir) {
    dup_biobank <- biobank0
    dup_biobank$Tube_Specific_ID[2] <- dup_biobank$Tube_Specific_ID[1]
    paths <- copy_inputs(case_dir, samples_df = samples0, biobank_df = dup_biobank)
    create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 2,
      seed = 123,
      biobank = paths$biobank,
      quiet = TRUE
    )
  }
)

add_case(
  name = "partial_inventory_does_not_process_unmatched_or_short_samples",
  notes = "Only fully satisfiable samples should be marked processed and consume tubes",
  runner = function(case_dir) {
    custom <- make_partial_inventory_data()
    paths <- copy_inputs(case_dir, samples_df = custom$samples, biobank_df = custom$biobank)
    samples_before <- read_csv_safe(paths$samples)
    biobank_before <- read_csv_safe(paths$biobank)
    pdf_file <- file.path(case_dir, "partial_inventory.pdf")
    samples_out <- file.path(case_dir, "samples_UPDATED.csv")
    biobank_out <- file.path(case_dir, "biobank_UPDATED.csv")

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 4,
      seed = 124,
      biobank = paths$biobank,
      used_col = "Used",
      replicates_col = "Replicates",
      tubes_per_sample = 2,
      pdf_file = pdf_file,
      write_updates = TRUE,
      update_in_place = FALSE,
      samples_out = samples_out,
      biobank_out = biobank_out,
      quiet = TRUE
    )

    processed_ids <- c(custom$ids[2], custom$ids[4])
    held_ids <- c(custom$ids[1], custom$ids[3])
    picked_tubes <- c("T2", "T3", "T4", "T5")
    held_tubes <- c("T1", "T6")

    check(nrow(res$batch) == 4L, "Expected all four custom samples to be selected.")
    check(nrow(res$biobank_picklist) == 4L, "Expected four picked tubes for fully satisfiable samples.")
    check(setequal(unique(res$biobank_picklist$Unique_ID), processed_ids), "Picklist should only contain fully satisfiable sample IDs.")
    check(!is.null(res$unmatched_samples), "Expected unmatched_samples for the sample with no Proteomic tubes.")
    check(!is.null(res$insufficient_tubes), "Expected insufficient_tubes for the sample with too few tubes.")
    check(identical(res$unmatched_samples$Unique_ID, custom$ids[3]), "Expected %s to be unmatched.", custom$ids[3])
    check(identical(res$insufficient_tubes$Unique_ID, custom$ids[1]), "Expected %s to be marked insufficient.", custom$ids[1])

    check(
      all(res$updated_samples[res$updated_samples$Unique_ID %in% processed_ids, "Proteomic"] == "No"),
      "Fully matched samples should be marked as No."
    )
    check(
      all(res$updated_samples[res$updated_samples$Unique_ID %in% held_ids, "Proteomic"] == "Yes"),
      "Unmatched or insufficient samples should remain Yes."
    )

    check(
      all(res$updated_biobank[res$updated_biobank$Tube_Specific_ID %in% picked_tubes, "Used"] == "Yes"),
      "Picked tubes should be marked Used == Yes in returned updated_biobank."
    )
    check(
      all(res$updated_biobank[res$updated_biobank$Tube_Specific_ID %in% held_tubes, "Used"] == "No"),
      "Unused or irrelevant tubes should remain Used == No."
    )

    if (!is.na(res$pdf)) {
      check(file.exists(samples_out), "Expected samples output file when PDF renders successfully.")
      check(file.exists(biobank_out), "Expected biobank output file when PDF renders successfully.")
      check(
        all(res$updated_biobank[res$updated_biobank$Tube_Specific_ID %in% picked_tubes, "Batch_sample_processed"] == basename(res$pdf)),
        "Picked tubes should be stamped with the rendered PDF file name."
      )
    } else {
      samples_after <- read_csv_safe(paths$samples)
      biobank_after <- read_csv_safe(paths$biobank)

      check(length(res$written_files) == 0L, "Expected no written files when PDF render fails.")
      check(same_df(samples_before, samples_after), "Samples input should remain unchanged when PDF render fails.")
      check(same_df(biobank_before, biobank_after), "Biobank input should remain unchanged when PDF render fails.")
      check(
        all(res$updated_biobank[res$updated_biobank$Tube_Specific_ID %in% picked_tubes, "Batch_sample_processed"] %in% c("", NA)),
        "Returned updated_biobank should not stamp Batch_sample_processed when PDF render fails."
      )
    }

    list(
      n_batch = nrow(res$batch),
      n_picklist = nrow(res$biobank_picklist),
      n_unmatched = nrow(res$unmatched_samples),
      n_insufficient = nrow(res$insufficient_tubes),
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = file.exists(samples_out),
      wrote_biobank = file.exists(biobank_out)
    )
  }
)

add_case(
  name = "no_matching_biobank_rows_keep_samples_eligible",
  notes = "When the biobank has no rows for the requested omic, samples should remain eligible",
  runner = function(case_dir) {
    no_match_biobank <- biobank0
    no_match_biobank$Application <- "Metabolomic"
    paths <- copy_inputs(case_dir, samples_df = samples0, biobank_df = no_match_biobank)
    pdf_file <- file.path(case_dir, "no_match_biobank.pdf")

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 4,
      seed = 125,
      biobank = paths$biobank,
      pdf_file = pdf_file,
      write_updates = FALSE,
      quiet = TRUE
    )

    selected_ids <- res$batch$Unique_ID
    check(is.null(res$biobank_picklist), "Expected NULL picklist when no biobank rows match the omic.")
    check(
      all(res$updated_samples[res$updated_samples$Unique_ID %in% selected_ids, "Proteomic"] == "Yes"),
      "Selected IDs should remain eligible when no biobank rows match the omic."
    )
    check(!is.null(res$unmatched_samples), "Expected unmatched_samples when no biobank rows match the omic.")
    check(nrow(res$unmatched_samples) == length(selected_ids), "All selected IDs should be unmatched when no biobank rows match.")

    list(
      n_batch = nrow(res$batch),
      n_picklist = 0L,
      n_unmatched = nrow(res$unmatched_samples),
      n_insufficient = 0L,
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = FALSE,
      wrote_biobank = FALSE
    )
  }
)

add_case(
  name = "missing_optional_biobank_cols_allowed_by_default",
  notes = "Used and Replicates columns should only be required when explicitly supplied",
  runner = function(case_dir) {
    biobank_min <- biobank0
    biobank_min$Used <- NULL
    biobank_min$Replicates <- NULL
    paths <- copy_inputs(case_dir, samples_df = samples0, biobank_df = biobank_min)
    pdf_file <- file.path(case_dir, "missing_optional_cols.pdf")

    res <- create_batch(
      samples = paths$samples,
      id_col = "Unique_ID",
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 4,
      seed = 126,
      biobank = paths$biobank,
      pdf_file = pdf_file,
      write_updates = FALSE,
      quiet = TRUE
    )

    check(!is.null(res$biobank_picklist), "Expected a picklist even when Used/Replicates columns are absent by default.")
    check(!"Used" %in% names(res$updated_biobank), "Used column should remain absent when it was not supplied in the input.")
    check(!"Replicates" %in% names(res$biobank_picklist), "Replicates should not appear in the picklist unless replicates_col is supplied.")

    list(
      n_batch = nrow(res$batch),
      n_picklist = nrow(res$biobank_picklist),
      n_unmatched = nrow(res$unmatched_samples %||% data.frame()),
      n_insufficient = nrow(res$insufficient_tubes %||% data.frame()),
      pdf_rendered = !is.na(res$pdf),
      wrote_samples = FALSE,
      wrote_biobank = FALSE
    )
  }
)

run_one_case <- function(case_obj, i, n_total, output_dir) {
  start_time <- Sys.time()
  warnings_seen <- character()
  messages_seen <- character()
  stdout_seen <- character()
  err_msg <- NA_character_
  status <- "ok"
  metrics <- list()

  safe_name <- gsub("[^A-Za-z0-9_]+", "_", case_obj$name)
  case_dir <- file.path(output_dir, sprintf("case_%04d_%s", i, safe_name))
  if (dir.exists(case_dir)) {
    unlink(case_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("[%s] %d/%d  %s\n", format(Sys.time(), "%H:%M:%S"), i, n_total, case_obj$name))

  stdout_seen <- capture.output(
    withCallingHandlers(
      tryCatch({
        metrics <- case_obj$runner(case_dir)
      }, error = function(e) {
        status <<- "error"
        err_msg <<- conditionMessage(e)
      }),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages_seen <<- c(messages_seen, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (any(grepl("The directory '.+' does not exist\\.", warnings_seen))) {
    status <- "error"
    if (is.na(err_msg) || !nzchar(err_msg)) {
      err_msg <- "Unexpected PDF output directory warning."
    }
  }
  if (any(grepl("Undefined control sequence\\.", messages_seen, fixed = FALSE))) {
    status <- "error"
    if (is.na(err_msg) || !nzchar(err_msg)) {
      err_msg <- "Unexpected LaTeX undefined control sequence during PDF rendering."
    }
  }
  if (any(grepl("Misplaced \\\\noalign\\.", messages_seen, fixed = FALSE))) {
    status <- "error"
    if (is.na(err_msg) || !nzchar(err_msg)) {
      err_msg <- "Unexpected LaTeX table-formatting error during PDF rendering."
    }
  }
  if (any(grepl("Environment .+ undefined\\.", messages_seen))) {
    status <- "error"
    if (is.na(err_msg) || !nzchar(err_msg)) {
      err_msg <- "Unexpected missing LaTeX environment during PDF rendering."
    }
  }
  passed <- if (isTRUE(case_obj$expect_error)) identical(status, "error") else identical(status, "ok")

  log_file <- file.path(output_dir, sprintf("case_%04d_%s.log", i, safe_name))
  writeLines(
    c(
      sprintf("case_id: %d", i),
      sprintf("name: %s", case_obj$name),
      sprintf("expect_error: %s", case_obj$expect_error),
      sprintf("status: %s", status),
      sprintf("pass: %s", passed),
      sprintf("elapsed_sec: %.3f", elapsed),
      sprintf("notes: %s", case_obj$notes),
      sprintf("case_dir: %s", case_dir),
      "",
      "STDOUT:",
      if (length(stdout_seen) > 0L) stdout_seen else "<none>",
      "",
      "MESSAGES:",
      if (length(messages_seen) > 0L) messages_seen else "<none>",
      "",
      "WARNINGS:",
      if (length(warnings_seen) > 0L) warnings_seen else "<none>",
      "",
      "ERROR:",
      if (!is.na(err_msg)) err_msg else "<none>"
    ),
    con = log_file
  )

  data.frame(
    case_id = i,
    name = case_obj$name,
    expect_error = case_obj$expect_error,
    status = status,
    pass = passed,
    elapsed_sec = round(elapsed, 3),
    n_batch = metrics$n_batch %||% NA_integer_,
    n_picklist = metrics$n_picklist %||% NA_integer_,
    n_unmatched = metrics$n_unmatched %||% NA_integer_,
    n_insufficient = metrics$n_insufficient %||% NA_integer_,
    pdf_rendered = metrics$pdf_rendered %||% NA,
    wrote_samples = metrics$wrote_samples %||% NA,
    wrote_biobank = metrics$wrote_biobank %||% NA,
    n_warnings = length(warnings_seen),
    n_messages = length(messages_seen),
    n_stdout = length(stdout_seen),
    error_message = ifelse(is.na(err_msg), "", err_msg),
    first_warning = clean_text(warnings_seen, n = 1L),
    first_message = clean_text(messages_seen, n = 1L),
    first_stdout = clean_text(stdout_seen, n = 1L),
    log_file = log_file,
    notes = case_obj$notes,
    stringsAsFactors = FALSE
  )
}

results <- vector("list", length(cases))
for (i in seq_along(cases)) {
  results[[i]] <- tryCatch(
    run_one_case(cases[[i]], i = i, n_total = length(cases), output_dir = OUTPUT_DIR),
    error = function(e) {
      data.frame(
        case_id = i,
        name = cases[[i]]$name,
        expect_error = cases[[i]]$expect_error,
        status = "runner_error",
        pass = FALSE,
        elapsed_sec = NA_real_,
        n_batch = NA_integer_,
        n_picklist = NA_integer_,
        n_unmatched = NA_integer_,
        n_insufficient = NA_integer_,
        pdf_rendered = NA,
        wrote_samples = NA,
        wrote_biobank = NA,
        n_warnings = NA_integer_,
        n_messages = NA_integer_,
        n_stdout = NA_integer_,
        error_message = conditionMessage(e),
        first_warning = "",
        first_message = "",
        first_stdout = "",
        log_file = "",
        notes = "Test runner internal error",
        stringsAsFactors = FALSE
      )
    }
  )
}

results_df <- do.call(rbind, results)

results_file <- file.path(OUTPUT_DIR, "create_batch_test_results.csv")
failures_file <- file.path(OUTPUT_DIR, "create_batch_test_failures.csv")
summary_file <- file.path(OUTPUT_DIR, "create_batch_test_summary.txt")
session_file <- file.path(OUTPUT_DIR, "sessionInfo.txt")

write.csv(results_df, results_file, row.names = FALSE)
write.csv(results_df[!results_df$pass, , drop = FALSE], failures_file, row.names = FALSE)

summary_lines <- c(
  sprintf("Total cases: %d", nrow(results_df)),
  sprintf("Passed: %d", sum(results_df$pass, na.rm = TRUE)),
  sprintf("Failed: %d", sum(!results_df$pass, na.rm = TRUE)),
  sprintf("Errors: %d", sum(results_df$status == 'error', na.rm = TRUE)),
  sprintf("Runner errors: %d", sum(results_df$status == 'runner_error', na.rm = TRUE)),
  sprintf("Cases with warnings: %d", sum(results_df$n_warnings > 0, na.rm = TRUE)),
  sprintf("Cases with messages: %d", sum(results_df$n_messages > 0, na.rm = TRUE)),
  sprintf("Cases with stdout: %d", sum(results_df$n_stdout > 0, na.rm = TRUE)),
  "",
  sprintf("Results CSV: %s", results_file),
  sprintf("Failures CSV: %s", failures_file),
  sprintf("Case logs directory: %s", OUTPUT_DIR)
)
writeLines(summary_lines, con = summary_file)
writeLines(capture.output(sessionInfo()), con = session_file)

cat("\n=== TEST SUMMARY ===\n")
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")

if (any(!results_df$pass)) {
  cat("\nTop failing cases:\n")
  failing <- results_df[!results_df$pass, c("case_id", "name", "status", "error_message"), drop = FALSE]
  print(utils::head(failing, 20))
}
