#' bisrDE: Bulk RNA-Seq Differential Expression Pipeline (BISR)
#'
#' DESeq2-based bulk RNA-seq differential expression pipeline for the
#' Bioinformatics Shared Resource at VCU Massey Comprehensive Cancer Center.
#'
#' Per-comparison DE, sample-level QC, GO + KEGG + Reactome + MSigDB Hallmark
#' enrichment, and a parameterised Quarto report. Designed to plug in
#' downstream of nf-core/rnaseq merged-counts output.
#'
#' @section Phase 5 in-progress notice:
#' The package skeleton is freshly scaffolded (Phase 5.1). Sub-phases 5.2-5.5
#' will migrate functions from the parent's `de.R` and `report_generator.R`
#' into per-domain files under `R/` (io, normalize, de, annotate, enrich,
#' plots_*, pipeline, cli). Phase 5.6 will move the report into
#' `inst/qmd/report.qmd` (Quarto). Until then, run-time pipeline still lives
#' in the parent's `de.R`.
#'
#' @keywords internal
#' @importFrom rlang .data :=
"_PACKAGE"

# Silence R CMD check NOTEs that come from S4 slot accessors and other
# bindings that lintr/check can't resolve through dispatch. These are NOT
# actually unbound — they're accessed via S4 `@` (e.g. `gse@result`) or
# present in the calling environment of dplyr/cli interpolations.
utils::globalVariables(c(
  "result",   # S4 slot of clusterProfiler GSEA result objects
  "PC1", "PC2", "PC3",  # plotly NSE inside `~PC1` formula refs
  "Group", "X", "Y"     # plotly + ggplot NSE refs
))
