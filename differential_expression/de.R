#!/usr/bin/env Rscript
#
# de.R — Bulk RNA-Seq Differential Expression Pipeline (thin wrapper).
#
# Phase 5.8 reduced this from a ~1880-line monolith to a ~70-line CLI
# driver. All pipeline logic now lives in the `bisrDE` package
# (./bisrDE/). This script just parses the same optparse flags as
# before and delegates to `bisrDE::run_pipeline()` + `bisrDE::generate_report()`.
#
# Usage (typically via run_analysis.sh which forwards args here):
#   Rscript de.R \
#     --counts <merged_counts.tsv> \
#     --samplesheet <samplesheet.csv> \
#     --outdir <out_dir> \
#     --runid <unique_id> \
#     --annotation <mouse|human> \
#     [--brs-ticket BRS-XXXX] \
#     [--id-type <ensembl|entrez|symbol>]

suppressPackageStartupMessages({
  library(optparse)
})

# ---------------------------------------------------------------------------
# Load the bisrDE package. Prefer the installed version (in container builds
# bisrDE is installed via remotes::install_local in dge_analysis.def). Fall
# back to devtools::load_all() against the sibling source dir for local
# development without an install step.
# ---------------------------------------------------------------------------
load_bisrDE <- function() {
  if (requireNamespace("bisrDE", quietly = TRUE)) {
    library(bisrDE)
    return(invisible(NULL))
  }
  # Find the directory containing THIS script (de.R) so the dev-mode
  # fallback works when invoked from a different cwd (e.g. a Nextflow
  # work dir several levels deep).
  script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[[1]]),
                            mustWork = FALSE))
    } else NULL
  }, error = function(e) NULL)

  candidates <- c(
    "bisrDE",
    "./bisrDE",
    file.path(getwd(), "bisrDE"),
    if (!is.null(script_dir)) file.path(script_dir, "bisrDE") else NULL
  )
  candidates <- Filter(Negate(is.null), candidates)

  for (p in candidates) {
    if (file.exists(file.path(p, "DESCRIPTION"))) {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        stop("bisrDE not installed and devtools not available for dev fallback. ",
             "Install bisrDE with: Rscript -e 'remotes::install_local(\"bisrDE\")'")
      }
      message(sprintf("[de.R] Dev mode: devtools::load_all('%s')", p))
      devtools::load_all(p, quiet = TRUE)
      return(invisible(NULL))
    }
  }
  stop("bisrDE not installed and source dir not found. ",
       "Install with: Rscript -e 'remotes::install_local(\"bisrDE\")'")
}
load_bisrDE()

# ---------------------------------------------------------------------------
# CLI flags (preserved verbatim from the pre-5.8 monolith).
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("-c", "--counts"),
    type = "character", default = NULL,
    help = "Required. Path to merged counts.tsv (typically rsem.merged.gene_counts.tsv from nf-core/rnaseq)."
  ),
  make_option(c("-s", "--samplesheet"),
    type = "character", default = NULL,
    help = "Required. Path to comma-delimited samplesheet.csv."
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "./output",
    help = "Output directory [default %default]."
  ),
  make_option(c("-r", "--runid"),
    type = "character", default = NULL,
    help = "Required. Unique name for this analysis run.",
    metavar = "character"
  ),
  make_option(c("-a", "--annotation"),
    type = "character", default = "mouse",
    help = "Genome annotation: 'mouse' or 'human' [default %default]."
  ),
  make_option(c("-b", "--brs-ticket"),
    type = "character", default = NULL,
    help = "Optional BRS ticket identifier (e.g. BRS-1234), rendered as a subtitle in the report."
  ),
  make_option(c("-i", "--id-type"),
    type = "character", default = "ensembl",
    help = "Gene identifier type used in the counts matrix rownames: 'ensembl' (default), 'entrez', or 'symbol'."
  ),
  make_option("--exclude-samples",
    type = "character", default = NULL,
    help = "Optional comma-separated SampleIDs to drop from the entire analysis (e.g. 'SRR1,SRR2')."
  ),
  make_option("--exclude-groups",
    type = "character", default = NULL,
    help = "Optional comma-separated GroupIDs to drop from the entire analysis."
  ),
  make_option("--include-contrasts",
    type = "character", default = NULL,
    help = "Optional comma-separated contrast columns to process exclusively (allowlist)."
  ),
  make_option("--exclude-contrasts",
    type = "character", default = NULL,
    help = "Optional comma-separated contrast columns to skip (denylist)."
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$runid) | is.null(opt$counts) |
    is.null(opt$samplesheet) | is.null(opt$annotation)) {
  print_help(opt_parser)
  stop("All required arguments must be supplied (input file)", call. = FALSE)
}

# optparse stores `--brs-ticket` and `--id-type` under their literal long-
# option names; coerce missing/NULL values to defaults so downstream code
# can treat them as simple strings.
brs_ticket <- if (is.null(opt[["brs-ticket"]])) "" else opt[["brs-ticket"]]
id_type    <- if (is.null(opt[["id-type"]]))    "ensembl" else opt[["id-type"]]

# Comma-separated list flags -> character vectors (NULL when unset/empty).
split_csv <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(NULL)
  vals <- trimws(strsplit(x, ",")[[1]])
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0) NULL else vals
}
exclude_samples   <- split_csv(opt[["exclude-samples"]])
exclude_groups    <- split_csv(opt[["exclude-groups"]])
include_contrasts <- split_csv(opt[["include-contrasts"]])
exclude_contrasts <- split_csv(opt[["exclude-contrasts"]])

# ---------------------------------------------------------------------------
# Run the pipeline + render the report.
# ---------------------------------------------------------------------------
res <- bisrDE::run_pipeline(
  counts_path       = opt$counts,
  samplesheet_path  = opt$samplesheet,
  outdir            = opt$outdir,
  runid             = opt$runid,
  annotation        = opt$annotation,
  id_type           = id_type,
  brs_ticket        = brs_ticket,
  exclude_samples   = exclude_samples,
  exclude_groups    = exclude_groups,
  include_contrasts = include_contrasts,
  exclude_contrasts = exclude_contrasts
)

bisrDE::generate_report(
  rds_path   = res$rds_path,
  output_dir = opt$outdir,
  brs_ticket = brs_ticket
)
