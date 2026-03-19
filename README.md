# batchcreator

`batchcreator` is a lightweight R package for creating balanced wet-lab sample
batches, matching selected samples to biobank inventory, and generating PDF
pick lists for laboratory workflows.

## What it includes

- Balanced batch selection with proportional or stratified sampling
- Optional biobank matching and tube pick-list generation
- Optional updates to sample and inventory status tables
- PDF batch reports with composition summaries, selected samples, and pick lists
- Example input files in `inst/extdata/`

## Install from GitHub

```r
install.packages("remotes")
remotes::install_github("bradleyalexward/batchcreator")
```

## Example

```r
library(batchcreator)

samples_path <- system.file("extdata", "samples.csv", package = "batchcreator")
biobank_path <- system.file("extdata", "biobank.csv", package = "batchcreator")

out <- create_batch(
  samples = samples_path,
  strata_cols = c("Sample_group", "Sex"),
  omic = "Proteomic",
  batch_size = 10,
  balance_strategy = "stratified",
  biobank = biobank_path,
  used_col = "Used",
  replicates_col = "Replicates",
  pdf_file = tempfile(fileext = ".pdf")
)
```

## Development notes

- Update the author, email, and license metadata in `DESCRIPTION` before
  publishing.
- Run `devtools::document()` to regenerate `man/` and `NAMESPACE` after future
  edits.
- Run `devtools::check()` before pushing a release tag.
