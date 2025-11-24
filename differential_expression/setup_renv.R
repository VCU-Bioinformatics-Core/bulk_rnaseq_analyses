#!/usr/bin/env Rscript
#
# setup_renv.R
# Bootstrap script to initialize renv and restore packages from dge_renv.lock
#
# Usage:
#   Rscript setup_renv.R
#   OR
#   source("setup_renv.R")
#

cat(strrep("=", 70), "\n", sep = "")
cat("Differential Expression Analysis - renv Setup\n")
cat(strrep("=", 70), "\n\n", sep = "")

# Check if renv is installed
if (!requireNamespace("renv", quietly = TRUE)) {
    cat("Installing renv package...\n")
    install.packages("renv", repos = "https://cran.r-project.org")
}

# Load renv
library(renv)

# Get the directory where this script is located
script_dir <- tryCatch(
    {
        # Try to get the script path when run with Rscript
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
            dirname(sub(needle, "", cmdArgs[match]))
        } else {
            # If run interactively or via source(), use current directory
            getwd()
        }
    },
    error = function(e) {
        getwd()
    }
)

setwd(script_dir)
cat("Working directory:", getwd(), "\n\n")

# Check if renv is already initialized
if (dir.exists("renv")) {
    cat("✓ renv already initialized\n")
    cat("  Re-activating project...\n\n")
    renv::activate()
} else {
    cat("Initializing renv project...\n")
    cat("This may take a few minutes on first run.\n\n")

    # Initialize renv with bare = TRUE to not snapshot immediately
    renv::init(bare = TRUE, restart = FALSE)
    cat("✓ renv initialized\n\n")
}

# Check if dge_renv.lock exists
if (file.exists("dge_renv.lock")) {
    cat("Found dge_renv.lock file\n")
    cat("Restoring packages from lockfile...\n")
    cat("NOTE: This will take 30-60 minutes on first run.\n")
    cat("      Required packages:\n")
    cat("      - Bioconductor packages (DESeq2, clusterProfiler, etc.)\n")
    cat("      - CRAN packages (tidyverse, plotly, etc.)\n\n")

    # Copy dge_renv.lock to renv.lock
    file.copy("dge_renv.lock", "renv.lock", overwrite = TRUE)

    # Restore packages from the lockfile
    cat("Starting package restoration...\n")
    renv::restore(prompt = FALSE)

    cat("\n")
    cat(strrep("=", 70), "\n", sep = "")
    cat("✓ Setup complete!\n")
    cat(strrep("=", 70), "\n\n", sep = "")
    cat("Packages have been installed and renv is configured.\n")
    cat("The .Rprofile will automatically activate renv when you start R.\n\n")
    cat("To verify installation:\n")
    cat("  1. Restart R\n")
    cat("  2. Check that renv activates automatically\n")
    cat("  3. Run: library(DESeq2)\n\n")
} else {
    cat("⚠ WARNING: dge_renv.lock not found!\n")
    cat("  Expected location:", file.path(getwd(), "dge_renv.lock"), "\n")
    cat("  renv is initialized but no packages have been installed.\n")
}

cat("Project library location:", renv::paths$library(), "\n")
cat("Cache location:", renv::paths$cache(), "\n\n")
