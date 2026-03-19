testthat::test_that("create_batch validates unique IDs", {
  df <- data.frame(
    Unique_ID = c("A", "A"),
    Sample_group = c("g1", "g1"),
    Sex = c("M", "F"),
    Proteomic = c("Yes", "Yes"),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    create_batch(
      samples = df,
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 1,
      biobank = NULL,
      quiet = TRUE
    ),
    "not unique"
  )
})

testthat::test_that("create_batch rejects fractional batch sizes", {
  df <- data.frame(
    Unique_ID = c("A", "B"),
    Sample_group = c("g1", "g1"),
    Sex = c("M", "F"),
    Proteomic = c("Yes", "Yes"),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    create_batch(
      samples = df,
      strata_cols = c("Sample_group", "Sex"),
      omic = "Proteomic",
      batch_size = 0.5,
      biobank = NULL,
      quiet = TRUE
    ),
    "positive whole number"
  )
})
