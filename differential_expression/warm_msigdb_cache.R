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

# msigdbr 25.x keeps one .rds per species in that cache and fetches each from
# Zenodo with a fixed 300 second limit and no retry (load_msigdb_df in the
# package's R/utils.R). Zenodo has served the 34 MB human file at 75 KB/s,
# which cannot finish in 300 seconds. So the files are fetched here first,
# with a longer limit and three attempts, into the path and name msigdbr
# expects; msigdbr then finds them cached and skips its own download. The URLs
# and checksums are read from the installed package's function body, so a
# version bump never leaves a stale copy here. If that body has another shape
# (a different msigdbr layout), nothing is prefetched and msigdbr downloads
# as before.
msigdb_sources <- function() {
  fn <- get0("load_msigdb_df", envir = asNamespace("msigdbr"), inherits = FALSE)
  if (!is.function(fn)) return(NULL)
  body <- deparse(fn)
  pick <- function(pattern) {
    hits <- grep(pattern, body, value = TRUE)
    sub(pattern, "\\1", hits)
  }
  urls <- pick('^\\s*file_url <- "([^"]+)"$')
  md5s <- pick('^\\s*file_md5 <- "([0-9a-f]{32})"$')
  if (length(urls) == 0 || length(urls) != length(md5s)) return(NULL)
  data.frame(url = urls, md5 = md5s, stringsAsFactors = FALSE)
}

prefetch <- function(src) {
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(nrow(src))) {
    dest <- file.path(cache, sub("\\?.*$", "", basename(src$url[i])))
    if (file.exists(dest) && unname(tools::md5sum(dest)) == src$md5[i]) {
      cat("  ", basename(dest), " already cached\n", sep = "")
      next
    }
    for (attempt in 1:3) {
      started <- Sys.time()
      good <- tryCatch({
        curl::curl_download(src$url[i], dest, quiet = TRUE,
                            handle = curl::new_handle(timeout = 3600, connecttimeout = 30))
        unname(tools::md5sum(dest)) == src$md5[i]
      }, error = function(e) {
        cat("  ", basename(dest), " attempt ", attempt, " failed: ", conditionMessage(e), "\n", sep = "")
        FALSE
      })
      if (good) {
        secs <- round(as.numeric(difftime(Sys.time(), started, units = "secs")))
        cat("  ", basename(dest), " downloaded (", round(file.size(dest) / 2^20), " MB, ", secs, " s)\n", sep = "")
        break
      }
      if (file.exists(dest)) file.remove(dest)
      if (attempt < 3) Sys.sleep(30)
    }
  }
}

src <- msigdb_sources()
if (is.null(src)) {
  cat("prefetch: msigdbr download sources not recognised; msigdbr downloads on its own\n\n")
} else {
  cat("prefetch from the package's own sources (", nrow(src), " files):\n", sep = "")
  prefetch(src)
  cat("\n")
}

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
