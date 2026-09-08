#!/usr/bin/env Rscript
# setup.R -- installs every package the notebooks in this repo need.
#
# Run once, from the repo root:
#   Rscript setup.R
#
# Needs a normal internet connection (CRAN + Bioconductor + GitHub).
# Written for R >= 4.3. Takes a while the first time -- Seurat and its
# dependency tree are large, and spacexr compiles C++ code.

options(warn = 1)

cran_pkgs <- c(
  "Seurat",        # v5, spatial analysis core
  "SeuratObject",
  "tidyverse",     # dplyr/ggplot2/tidyr/purrr/readr/stringr/forcats/tibble
  "patchwork",     # combining ggplot panels
  "hdf5r",         # reads the .h5 count matrices
  "data.table",    # fast I/O for the GSEA results table
  "pheatmap",
  "arrow",         # reads Visium HD's tissue_positions.parquet
  "qs2",           # fast object serialization (saving the RCTD reference etc.)
  "msigdbr",       # MSigDB gene sets for GSEA, as an R data frame
  "scatterpie",    # per-spot cell-type composition pie charts
  "BiocManager",
  "remotes"
)

bioc_pkgs <- c(
  "fgsea",         # gene set enrichment
  "SPOTlight",     # only used here for its plotSpatialScatterpie() helper
  "glmGamPoi"      # optional, but SCTransform runs noticeably faster with it --
                    # worth having given how many times this pipeline calls SCTransform
)

install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    cat("All of:", paste(pkgs, collapse = ", "), "-- already installed.\n")
    return(invisible())
  }
  cat("Installing:", paste(missing, collapse = ", "), "\n")
  installer(missing)
}

install_if_missing(cran_pkgs, function(p) install.packages(p, Ncpus = max(1, parallel::detectCores() - 1)))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_if_missing(bioc_pkgs, function(p) BiocManager::install(p, update = FALSE, ask = FALSE))

# spacexr (RCTD) isn't on CRAN or Bioconductor -- GitHub only.
if (!requireNamespace("spacexr", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  cat("Installing spacexr (RCTD) from GitHub -- this compiles C++ code, give it a few minutes...\n")
  remotes::install_github("dmcable/spacexr", upgrade = "never")
}

# presto (GitHub only) makes FindMarkers/FindAllMarkers's Wilcoxon test
# much faster -- optional, but this pipeline calls those functions on
# 100K+ bins, where the difference is real. Seurat auto-detects and uses
# it if present; nothing else to configure.
if (!requireNamespace("presto", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  cat("Installing presto (optional, speeds up FindMarkers) from GitHub...\n")
  tryCatch(
    remotes::install_github("immunogenomics/presto", upgrade = "never"),
    error = function(e) message("presto install failed (non-fatal, just slower): ", conditionMessage(e))
  )
}

cat("\n== Verifying ==\n")
required_pkgs <- c(cran_pkgs, bioc_pkgs, "spacexr")
optional_pkgs <- c("glmGamPoi", "presto")  # already covered above, called out separately below

req_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = required_pkgs, installed = req_ok), row.names = FALSE)

opt_ok <- vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)
cat("\nOptional speed-ups:\n")
print(data.frame(package = optional_pkgs, installed = opt_ok), row.names = FALSE)

if (all(req_ok)) {
  cat("\nAll required packages installed. Next: bash download_data.sh\n")
  if (!all(opt_ok)) cat("(Optional speed-up packages above are missing -- fine to skip, just slower.)\n")
} else {
  cat("\nSome REQUIRED packages failed to install -- see messages above. Re-run this",
      "script after resolving them (already-installed packages are skipped).\n")
}
