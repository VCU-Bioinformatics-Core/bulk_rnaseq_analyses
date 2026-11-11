#!/usr/bin/env Rscript
#
# synthetic_counts.R: write a seeded synthetic bulk RNA-seq data set with a
# known set of differentially expressed genes.
#
# The data set is the shared input for validating one runtime against another
# (container versus conda, conda versus renv): generate it once, run the same
# files through each runtime, then compare the runs with compare_runs.R. The
# truth file lets a reader check that a run recovers the planted signal.
#
#   Rscript tests/synthetic_counts.R --out <dir>
#   Rscript tests/synthetic_counts.R --out <dir> --genes 5000 --per-group 4 \
#                                    --de-fraction 0.05 --seed 1
#
# Options:
#   --out <dir>          output directory, created if missing         (required)
#   --genes <n>          number of genes                              [20000]
#   --per-group <n>      samples per group                            [3]
#   --de-fraction <x>    fraction of genes that are truly DE          [0.1]
#   --seed <n>           seed for the single set.seed() call          [20261015]
#
# Files written to --out:
#   counts.tsv       gene_id, then one integer count column per sample
#   samplesheet.csv  SampleID, GroupID, treated_vs_control (1 treated, 0 control)
#   truth.csv        gene_id, is_de, log2fc_true (treated over control)
#   manifest.txt     seed, sizes, org.Hs.eg.db and R versions
#
# Gene IDs are human Ensembl IDs drawn from the org.Hs.eg.db keys that map to
# both an ENTREZID and a SYMBOL, so annotation and every enrichment backend have
# real genes to work with. The ID universe depends on the org.Hs.eg.db version:
# generate once and copy the written files to every environment under test
# instead of regenerating in each. With the same arguments and the same
# org.Hs.eg.db version, two runs write byte-identical counts.tsv and truth.csv.
#
# Needs AnnotationDbi and org.Hs.eg.db (the renv library or the conda
# environment provides both).

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}
has_flag <- function(flag) flag %in% args

if (has_flag("-h") || has_flag("--help") || length(args) == 0) {
  cat("
synthetic_counts.R: seeded synthetic counts with a known DE gene set

Required:
  --out <dir>          output directory

Optional:
  --genes <n>          number of genes                      [20000]
  --per-group <n>      samples per group                    [3]
  --de-fraction <x>    fraction of truly DE genes           [0.1]
  --seed <n>           random seed                          [20261015]

Writes counts.tsv, samplesheet.csv, truth.csv and manifest.txt to --out.
")
  quit(status = 0)
}

out_dir <- get_arg("--out")
if (is.null(out_dir)) {
  message("ERROR: --out is required (see --help)")
  quit(status = 3)
}
n_genes   <- as.integer(get_arg("--genes", "20000"))
per_group <- as.integer(get_arg("--per-group", "3"))
de_frac   <- as.numeric(get_arg("--de-fraction", "0.1"))
seed      <- as.integer(get_arg("--seed", "20261015"))
if (anyNA(c(n_genes, per_group, de_frac, seed)) || n_genes < 1L ||
    per_group < 2L || de_frac < 0 || de_frac > 1) {
  message("ERROR: --genes >= 1, --per-group >= 2 and 0 <= --de-fraction <= 1 are required")
  quit(status = 3)
}

for (pkg in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("ERROR: package '", pkg, "' is not installed. Run this script in the ",
            "renv library or the conda environment of the pipeline.")
    quit(status = 3)
  }
}

# ---------------------------------------------------------------------------
# Gene universe
# ---------------------------------------------------------------------------
db  <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
ens <- sort(AnnotationDbi::keys(db, keytype = "ENSEMBL"))
map <- suppressMessages(AnnotationDbi::select(
  db, keys = ens, keytype = "ENSEMBL", columns = c("ENTREZID", "SYMBOL")))
map <- map[!is.na(map$ENTREZID) & !is.na(map$SYMBOL), , drop = FALSE]
universe <- sort(unique(map$ENSEMBL))
if (length(universe) < n_genes) {
  message("ERROR: --genes ", n_genes, " exceeds the ", length(universe),
          " Ensembl IDs with an ENTREZID and a SYMBOL in org.Hs.eg.db")
  quit(status = 3)
}

# ---------------------------------------------------------------------------
# Simulation: one set.seed(), then every random draw in a fixed order
# ---------------------------------------------------------------------------
set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")

genes <- sort(sample(universe, n_genes))

group   <- rep(c("control", "treated"), each = per_group)
samples <- c(paste0("control_", seq_len(per_group)),
             paste0("treated_", seq_len(per_group)))
treated <- as.integer(group == "treated")

base_mean <- pmax(1, rlnorm(n_genes, meanlog = 4, sdlog = 1.5))
size_fac  <- runif(length(samples), min = 0.7, max = 1.4)
disp      <- 0.05 + 2 / base_mean

n_de   <- as.integer(round(de_frac * n_genes))
is_de  <- rep(FALSE, n_genes)
is_de[sample.int(n_genes, n_de)] <- TRUE
lfc    <- numeric(n_genes)
lfc[is_de] <- sample(c(-1, 1), n_de, replace = TRUE) *
  runif(n_de, min = 1, max = 3)

mu <- outer(base_mean, size_fac) * 2^outer(lfc, treated)
counts <- matrix(
  stats::rnbinom(length(mu), mu = as.vector(mu), size = rep(1 / disp, times = length(samples))),
  nrow = n_genes, dimnames = list(genes, samples))
storage.mode(counts) <- "integer"

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_df <- data.frame(gene_id = genes, counts, check.names = FALSE,
                        stringsAsFactors = FALSE)
utils::write.table(counts_df, file.path(out_dir, "counts.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)

sheet <- data.frame(SampleID = samples, GroupID = group,
                    treated_vs_control = ifelse(treated == 1L, 1L, 0L),
                    stringsAsFactors = FALSE)
utils::write.csv(sheet, file.path(out_dir, "samplesheet.csv"),
                 quote = FALSE, row.names = FALSE)

truth <- data.frame(gene_id = genes, is_de = is_de,
                    log2fc_true = signif(lfc, 6), stringsAsFactors = FALSE)
utils::write.csv(truth, file.path(out_dir, "truth.csv"),
                 quote = FALSE, row.names = FALSE)

writeLines(c(
  paste0("seed: ", seed),
  paste0("genes: ", n_genes),
  paste0("samples_per_group: ", per_group),
  paste0("de_fraction: ", de_frac),
  paste0("de_genes: ", n_de),
  paste0("org.Hs.eg.db: ", as.character(utils::packageVersion("org.Hs.eg.db"))),
  paste0("R: ", R.version.string)
), file.path(out_dir, "manifest.txt"))

cat("wrote ", n_genes, " genes x ", length(samples), " samples (",
    n_de, " DE) to ", out_dir, "\n", sep = "")
