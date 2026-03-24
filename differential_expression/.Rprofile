# .Rprofile for differential expression analysis
# This file automatically activates renv when R starts in this directory

# Source renv/activate.R if it exists
if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
    cat("✓ renv activated for differential expression analysis\n")
    cat("  Using R", as.character(getRversion()), "\n")
    cat("  Project library:", renv::paths$library(), "\n\n")
} else {
    cat("⚠ renv not initialized. Run: source('setup_renv.R')\n")
}
