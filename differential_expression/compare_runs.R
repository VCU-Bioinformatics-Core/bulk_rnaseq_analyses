#!/usr/bin/env Rscript
#
# compare_runs.R — diff two bisrDE runs to decide whether an environment change
# is safe.
#
# The intended use is validating the conda/bioconda environment (Bioconductor
# 3.18) against the renv reference (Bioconductor 3.16): run the same data
# through both, then diff the DE results before trusting the new environment
# for production analyses.
#
#   Rscript compare_runs.R --a <reference_outdir> --b <candidate_outdir>
#   Rscript compare_runs.R --a old/ --b new/ --padj 0.05 --fold-change 1.5 \
#                          --out validation.md
#
# Exit status: 0 = PASS, 1 = WARN, 2 = FAIL, 3 = usage/IO error. So it can gate
# a deployment in CI.
#
# Deliberately depends on BASE R ONLY — no bisrDE, no tidyverse — so it runs in
# either environment (or a bare R) without becoming part of what it validates.

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}
has_flag <- function(flag) flag %in% args

if (has_flag("-h") || has_flag("--help") || length(args) == 0) {
  cat("
compare_runs.R — diff two bisrDE runs (e.g. renv reference vs conda candidate)

Required:
  --a <dir>            reference run output directory
  --b <dir>            candidate run output directory

Optional:
  --padj <num>         significance cutoff for DEG sets      [0.05]
  --fold-change <num>  linear fold-change cutoff for DEG sets [1.5]
  --out <file>         also write a markdown report
  --min-cor <num>      min log2FC Spearman correlation to PASS [0.99]
  --min-jaccard <num>  min DEG-set Jaccard to PASS             [0.95]

Exit: 0 PASS, 1 WARN, 2 FAIL, 3 usage/IO error
")
  quit(status = 0)
}

dir_a <- get_arg("--a")
dir_b <- get_arg("--b")
if (is.null(dir_a) || is.null(dir_b)) {
  message("ERROR: --a and --b are both required (see --help)")
  quit(status = 3)
}
padj_cut    <- as.numeric(get_arg("--padj", "0.05"))
fc_cut      <- as.numeric(get_arg("--fold-change", "1.5"))
lfc_cut     <- log2(fc_cut)
out_file    <- get_arg("--out", NA_character_)
min_cor     <- as.numeric(get_arg("--min-cor", "0.99"))
min_jaccard <- as.numeric(get_arg("--min-jaccard", "0.95"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
de_dir <- function(root) file.path(root, "data", "de_data")

# Map comparison name -> DE csv path.
de_files <- function(root) {
  d <- de_dir(root)
  if (!dir.exists(d)) return(setNames(character(0), character(0)))
  f <- list.files(d, pattern = "^DESeq2_.*\\.csv$", full.names = TRUE)
  setNames(f, sub("\\.csv$", "", sub("^DESeq2_", "", basename(f))))
}

# The gene key is whichever ID column both tables share. Prefer the input-ID
# column, since that is what the counts matrix was keyed by.
pick_key <- function(x, y) {
  for (k in c("ENSEMBL_ID", "SYMBOL_ID", "ENTREZID_ID", "SYMBOL", "ENTREZID")) {
    if (k %in% names(x) && k %in% names(y)) return(k)
  }
  NULL
}

deg_set <- function(df, key) {
  ok <- !is.na(df$padj) & df$padj < padj_cut & abs(df$log2FoldChange) >= lfc_cut
  unique(as.character(df[[key]][ok]))
}

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

fmt <- function(x, d = 4) if (is.na(x)) "n/a" else formatC(x, format = "f", digits = d)

# Enrichment tables, keyed by backend directory.
enrich_terms <- function(root, comparison) {
  backends <- c(GO = "gsea_data", KEGG = "kegg_data",
                Reactome = "reactome_data", Hallmark = "hallmark_data")
  out <- list()
  for (nm in names(backends)) {
    d <- file.path(root, "data", backends[[nm]])
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = paste0(comparison, "\\.csv$"), full.names = TRUE)
    if (length(f) == 0) next
    tab <- tryCatch(utils::read.csv(f[1], stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(tab)) next
    idcol <- intersect(c("ID", "Description"), names(tab))
    if (length(idcol) == 0) next
    out[[nm]] <- unique(as.character(tab[[idcol[1]]]))
  }
  out
}

# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------
fa <- de_files(dir_a)
fb <- de_files(dir_b)

if (length(fa) == 0 || length(fb) == 0) {
  message("ERROR: no DESeq2_*.csv found under ", de_dir(dir_a),
          " and/or ", de_dir(dir_b))
  quit(status = 3)
}

only_a <- setdiff(names(fa), names(fb))
only_b <- setdiff(names(fb), names(fa))
shared <- intersect(names(fa), names(fb))

lines   <- character(0)
say <- function(...) {
  s <- paste0(...)
  cat(s, "\n", sep = "")
  lines <<- c(lines, s)
}

say("# bisrDE run comparison")
say("")
say("- reference (A): `", dir_a, "`")
say("- candidate (B): `", dir_b, "`")
say("- thresholds: padj < ", padj_cut, ", |log2FC| >= ", fmt(lfc_cut, 3),
    " (", fc_cut, "-fold)")
say("- PASS gates: log2FC Spearman >= ", min_cor, ", DEG Jaccard >= ", min_jaccard)
say("")

status <- 0L   # 0 pass, 1 warn, 2 fail
bump <- function(s) status <<- max(status, s)

if (length(only_a) || length(only_b)) {
  bump(2L)
  say("## Comparison sets differ")
  if (length(only_a)) say("- **only in A**: ", paste(only_a, collapse = ", "))
  if (length(only_b)) say("- **only in B**: ", paste(only_b, collapse = ", "))
  say("")
}

for (cmp in shared) {
  A <- tryCatch(utils::read.csv(fa[[cmp]], stringsAsFactors = FALSE), error = function(e) NULL)
  B <- tryCatch(utils::read.csv(fb[[cmp]], stringsAsFactors = FALSE), error = function(e) NULL)
  say("## ", cmp)
  if (is.null(A) || is.null(B)) {
    bump(2L); say("- **FAIL** — could not read one of the tables"); say(""); next
  }
  key <- pick_key(A, B)
  if (is.null(key)) {
    bump(2L); say("- **FAIL** — no shared gene-ID column"); say(""); next
  }

  A <- A[!duplicated(A[[key]]), , drop = FALSE]
  B <- B[!duplicated(B[[key]]), , drop = FALSE]
  common <- intersect(as.character(A[[key]]), as.character(B[[key]]))
  ia <- match(common, as.character(A[[key]]))
  ib <- match(common, as.character(B[[key]]))

  n_union <- length(union(as.character(A[[key]]), as.character(B[[key]])))
  gene_cov <- if (n_union == 0) NA_real_ else length(common) / n_union

  lfc_a <- suppressWarnings(as.numeric(A$log2FoldChange[ia]))
  lfc_b <- suppressWarnings(as.numeric(B$log2FoldChange[ib]))
  padj_a <- suppressWarnings(as.numeric(A$padj[ia]))
  padj_b <- suppressWarnings(as.numeric(B$padj[ib]))

  ok <- is.finite(lfc_a) & is.finite(lfc_b)
  rho <- if (sum(ok) > 2) suppressWarnings(stats::cor(lfc_a[ok], lfc_b[ok], method = "spearman")) else NA_real_
  pear <- if (sum(ok) > 2) suppressWarnings(stats::cor(lfc_a[ok], lfc_b[ok])) else NA_real_
  maxd <- if (any(ok)) max(abs(lfc_a[ok] - lfc_b[ok])) else NA_real_
  medd <- if (any(ok)) stats::median(abs(lfc_a[ok] - lfc_b[ok])) else NA_real_

  okp <- is.finite(padj_a) & is.finite(padj_b)
  prho <- if (sum(okp) > 2) suppressWarnings(stats::cor(padj_a[okp], padj_b[okp], method = "spearman")) else NA_real_

  da <- deg_set(A, key); db <- deg_set(B, key)
  jac <- jaccard(da, db)
  lost   <- setdiff(da, db)
  gained <- setdiff(db, da)

  say("- genes: A=", nrow(A), " B=", nrow(B), " shared=", length(common),
      " (", fmt(100 * gene_cov, 2), "% of union)")
  say("- log2FC: Spearman=", fmt(rho), " Pearson=", fmt(pear),
      " median|Δ|=", fmt(medd), " max|Δ|=", fmt(maxd))
  say("- padj: Spearman=", fmt(prho))
  say("- DEGs: A=", length(da), " B=", length(db),
      " shared=", length(intersect(da, db)), " Jaccard=", fmt(jac))
  if (length(lost))   say("  - lost in B (", length(lost), "): ",
                          paste(utils::head(lost, 10), collapse = ", "),
                          if (length(lost) > 10) " …" else "")
  if (length(gained)) say("  - gained in B (", length(gained), "): ",
                          paste(utils::head(gained, 10), collapse = ", "),
                          if (length(gained) > 10) " …" else "")

  ea <- enrich_terms(dir_a, cmp); eb <- enrich_terms(dir_b, cmp)
  for (bk in union(names(ea), names(eb))) {
    ta <- if (is.null(ea[[bk]])) character(0) else ea[[bk]]
    tb <- if (is.null(eb[[bk]])) character(0) else eb[[bk]]
    say("- ", bk, " terms: A=", length(ta), " B=", length(tb),
        " Jaccard=", fmt(jaccard(ta, tb)))
  }

  verdict <- "PASS"
  if (!is.na(rho) && rho < min_cor)      { verdict <- "FAIL"; bump(2L) }
  if (jac < min_jaccard)                 { if (verdict == "PASS") verdict <- "WARN"; bump(1L) }
  if (!is.na(gene_cov) && gene_cov < 0.99) { verdict <- "FAIL"; bump(2L) }
  say("- **", verdict, "**")
  say("")
}

overall <- c("PASS", "WARN", "FAIL")[status + 1L]
say("## Overall: ", overall)
if (status == 0L) {
  say("Results are equivalent within tolerance — the candidate environment looks safe.")
} else if (status == 1L) {
  say("Differences are small but real. Inspect the gained/lost DEGs before deploying.")
} else {
  say("Results differ materially. Do NOT treat the two environments as interchangeable.")
}

if (!is.na(out_file)) {
  writeLines(lines, out_file)
  cat("\nwrote ", out_file, "\n", sep = "")
}

quit(status = status)
