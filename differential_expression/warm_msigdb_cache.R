#!/usr/bin/env Rscript
#
# warm_msigdb_cache.R — pre-download the MSigDB gene sets for offline use.
#
# msigdbr >= 24 fetches its gene-set archive from Zenodo the first time it is
# used. HPC compute nodes usually have no outbound internet, so without a warm
# cache every run fails with:
#
#     MSigDB Hallmark error: Timeout was reached [zenodo.org]
#
# ...and the report silently contains no Hallmark enrichment. The archive is
# cached on disk and the download is skipped entirely once it is present, so
# running this ONCE on a machine with internet (e.g. a login node) fixes it
# permanently.
#
#     Rscript warm_msigdb_cache.R
#
# Deliberately requires nothing but R — no `just`, no conda tooling — because
# the machines that need it are exactly the ones without a dev toolchain.
#
# If $HOME is not shared with the compute nodes, put the cache somewhere that
# is, and set the same variable in your job script:
#
#     export R_USER_CACHE_DIR=/lustre/home/<lab>/rcache
#     Rscript warm_msigdb_cache.R

# NOTE: the conda-vs-renv library decision happens in .Rprofile, which has
# already run by the time this script starts — setting the environment variable
# here would be too late to matter. .Rprofile skips renv when $CONDA_PREFIX is
# set, so an activated conda environment is handled automatically.

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("Package 'msigdbr' is not available to this R.\n",
       "  R:        ", R.home(), "\n",
       "  libPaths: ", paste(.libPaths(), collapse = ", "), "\n",
       "If you are using the conda environment, activate it first; if R still\n",
       "cannot see it, renv's autoloader is shadowing the environment — see\n",
       "the HPC section of README.md.",
       call. = FALSE)
}

cache <- tools::R_user_dir("msigdbr", "cache")
cat("cache directory: ", cache, "\n", sep = "")
cat("(override with the R_USER_CACHE_DIR environment variable)\n\n")

species <- c("Homo sapiens", "Mus musculus")
ok <- TRUE

for (sp in species) {
  res <- tryCatch({
    n <- nrow(msigdbr::msigdbr(species = sp, collection = "H"))
    sprintf("  %-14s OK  (%d rows in the H collection)", sp, n)
  }, error = function(e) {
    ok <<- FALSE
    sprintf("  %-14s FAILED: %s", sp, conditionMessage(e))
  })
  cat(res, "\n", sep = "")
}

cat("\ncached files:\n")
files <- list.files(cache)
if (length(files) == 0) {
  cat("  (none — the download did not complete)\n")
} else {
  cat(paste0("  ", files, collapse = "\n"), "\n", sep = "")
}

if (!ok) {
  cat("\nThe download failed. This machine needs outbound internet access;\n",
      "run it on a login node, or copy an already-warmed cache directory here.\n",
      sep = "")
  quit(status = 1)
}

cat("\nDone. Compute nodes will now read this cache instead of downloading.\n")
