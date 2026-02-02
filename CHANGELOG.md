# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2026-02-02

### Added

- **Self-contained R environment** using `renv` for reproducible package management
- **Singularity/Apptainer container support** for Linux server execution
- **Launcher script** (`run_analysis.sh`) that auto-detects platform and execution method
- **Automated HTML report generation** with interactive visualizations
- **PCA visualization** with 2D and 3D interactive plots (Plotly)
- **GSEA (Gene Set Enrichment Analysis)** with GO term enrichment and dotplots

### Fixed

- **Report generator duplicate chunk error** - Fixed `glue` vectorization issue when comparisons have multiple experimental groups
- **GSEA `process_gsea` function** - Removed duplicate `tryCatch` blocks and corrected variable references
- **Figure paths in HTML reports** - Changed from absolute to relative paths so images display correctly
- **RDS file path resolution** - Used `normalizePath()` to ensure correct path handling during report generation
- **Missing `glue` library** - Added required import in `report_generator.R`

### Changed

- Improved error handling throughout the analysis pipeline
- Enhanced logging and progress output during analysis
- Updated volcano plot and heatmap generation for better visualization

## [Unreleased]

### Planned

- Support for additional annotation databases
- Batch correction options
- Enhanced QC metrics reporting

---

## Version History

| Version | Date       | Description                          |
| ------- | ---------- | ------------------------------------ |
| 1.3.0   | 2026-02-02 | Stable release with full DE pipeline |
