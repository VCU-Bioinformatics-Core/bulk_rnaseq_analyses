# .Rprofile for differential expression analysis
# This file automatically activates renv when R starts in this directory

# Skip renv entirely when a conda environment owns the R library, or when the
# renv autoloader has been switched off explicitly. Without this, sourcing
# activate.R below repoints .libPaths() at the (empty) project library and
# hides everything conda installed — the symptom is a confusing
# "there is no package called '<x>'" for packages that are demonstrably there.
.bisr_skip_renv <- nzchar(Sys.getenv("CONDA_PREFIX")) ||
  identical(toupper(Sys.getenv("RENV_CONFIG_AUTOLOADER_ENABLED")), "FALSE")

if (.bisr_skip_renv) {
  cat("• renv skipped (conda environment in use)\n")
  rm(.bisr_skip_renv)
} else

# Source renv/activate.R if it exists
if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
    cat("✓ renv activated for differential expression analysis\n")
    cat("  Using R", as.character(getRversion()), "\n")
    cat("  Project library:", renv::paths$library(), "\n\n")
} else {
    cat("⚠ renv not initialized. Run: source('setup_renv.R')\n")
}
