#!/usr/bin/env Rscript
#
# update_lockfile.R
# Updates the lockfile to use current Bioconductor version
#

cat("======================================================================\n")
cat("Updating lockfile for current R/Bioconductor versions\n")
cat("======================================================================\n\n")

# Check R version
r_version <- getRversion()
cat("Current R version:", as.character(r_version), "\n")

# Determine appropriate Bioconductor version
if (r_version >= "4.4.0") {
    bioc_version <- "3.19"
} else if (r_version >= "4.3.0") {
    bioc_version <- "3.18"
} else if (r_version >= "4.2.0") {
    bioc_version <- "3.16"
} else {
    stop("R version too old. Please upgrade to R >= 4.2.0")
}

cat("Will use Bioconductor version:", bioc_version, "\n\n")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

library(BiocManager)

# Set Bioconductor version
BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE)

# Install renv if needed
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}

library(renv)

# Initialize renv if not already done
if (!dir.exists("renv")) {
    cat("Initializing renv...\n")
    renv::init(bare = TRUE, restart = FALSE)
}

cat("\nInstalling required packages...\n")
cat("This will take 30-45 minutes.\n\n")

# Define all required packages
packages <- c(
    # CRAN packages
    "here", "dplyr", "data.table", "tidyverse", "janitor", "scales",
    "ggrepel", "readr", "DT", "ggplot2", "gplots", "RColorBrewer",
    "purrr", "plotly", "optparse", "htmlwidgets", "rmarkdown", "knitr",

    # Bioconductor packages
    "BiocManager", "DESeq2", "apeglm", "edgeR", "AnnotationDbi",
    "clusterProfiler", "enrichplot", "org.Mm.eg.db", "org.Hs.eg.db"
)

# Install packages
for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        if (pkg %in% c(
            "DESeq2", "apeglm", "edgeR", "AnnotationDbi", "clusterProfiler",
            "enrichplot", "org.Mm.eg.db", "org.Hs.eg.db"
        )) {
            BiocManager::install(pkg, ask = FALSE, update = FALSE)
        } else {
            install.packages(pkg)
        }
    }
}

cat("\nCreating snapshot...\n")
renv::snapshot(prompt = FALSE)

cat("\n======================================================================\n")
cat("✓ Lockfile updated successfully!\n")
cat("======================================================================\n\n")
cat("The renv.lock file now uses Bioconductor", bioc_version, "\n")
cat("Compatible with R", as.character(r_version), "\n\n")
cat("Next steps:\n")
cat("  1. Test the installation: Rscript -e 'library(DESeq2)'\n")
cat("  2. Run analysis: bash run_analysis.sh --help\n\n")
