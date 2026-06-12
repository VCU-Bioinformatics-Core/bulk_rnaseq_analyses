# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.0] - 2026-06-12

User-experience release: selective sample/group/contrast inclusion, custom
plot labels, and (in progress) a charmbracelet-powered interactive CLI.

### Added

- **Sample / group exclusion** — drop samples from the *entire* analysis
  (contrasts, normalization, DESeq2, plots) via four complementary
  mechanisms, all unioned:
  - `--exclude-samples SRR1,SRR2` — drop by `SampleID` (ad-hoc, no
    samplesheet edit).
  - `--exclude-groups GroupA,GroupB` — drop every sample in those groups.
  - Optional `Exclude` samplesheet column — declarative per-sample drop
    (`1` / `TRUE` / `yes`), a permanent record kept with the data.
  - `filter_samplesheet()` (new exported fn) implements the union and
    errors if fewer than two samples survive.
- **Contrast selection** — process a subset of contrasts without editing
  the samplesheet:
  - `--include-contrasts a_vs_b` — allowlist.
  - `--exclude-contrasts c_vs_d` — denylist.
  - `parse_contrasts()` gained `include_contrasts` / `exclude_contrasts`
    args and is now metadata-aware (skips `SampleID` / `GroupID` /
    `Exclude` / `DisplayName` columns by name, so a metadata column between
    `GroupID` and the contrasts is never mistaken for a contrast).
- **Custom plot labels** — optional `DisplayName` samplesheet column. When
  present, all sample-labelled figures (correlation + vst-distance
  heatmaps, library-size barplot, hclust dendrogram, static + interactive
  PCA, per-comparison DE heatmaps) show the display label instead of the
  SRA accession / SampleID. The matrix join key (`make.names(SampleID)`)
  is untouched — only labels change. Blank `DisplayName` cells fall back to
  the SampleID. New internal `.display_lookup()` centralises the mapping.
- **Nextflow params** `--exclude_samples`, `--exclude_groups`,
  `--include_contrasts`, `--exclude_contrasts` forwarded through the
  `nf-module` wrapper.
- **Interactive launcher** `run_interactive.sh` (charmbracelet, Phase A) —
  a guided, styled front-door that collects every parameter and offers
  interactive group / sample / contrast multi-select (wiring into the new
  selection features), then hands off to `run_analysis.sh`. Uses
  [`gum`](https://github.com/charmbracelet/gum) (lipgloss borders + huh
  forms + bubbles spinner/multi-select), [`glow`](https://github.com/charmbracelet/glow)
  (terminal-markdown run summary) and [`freeze`](https://github.com/charmbracelet/freeze)
  (CLI screenshots for docs) when installed; degrades gracefully to plain
  prompts otherwise. Scriptable via `BISR_*` env vars + `--print-cmd`.
  A dedicated Go TUI (bubbletea) is scoped as an optional Phase B.
- New testthat files `test-exclude.R` (11 tests) and `test-display-names.R`
  (7 tests) covering exclusion union, contrast filtering, metadata-column
  skipping, and the display-label join-key invariant.

### Changed

- `run_pipeline()` gained `exclude_samples` / `exclude_groups` /
  `include_contrasts` / `exclude_contrasts` parameters (all `NULL` by
  default — behaviour unchanged when unset).

## [1.4.0] - 2026-05-01

The biggest release since the pipeline was first published. The monolithic
`de.R` (~1880 lines) and glue-built `report_generator.R` (~625 lines) have
been refactored into a real R package (`bisrDE`) with a thin CLI wrapper, a
Quarto-templated report, a Nextflow DSL2 module, and a testthat suite.

### Added

- **`bisrDE` R package** at `differential_expression/bisrDE/` with 33
  exported functions across 11 R files (io / normalize / de / annotate /
  enrich / 4× plots_* / pipeline / cli / utils / report). Installable via
  `remotes::install_local("bisrDE")`. Roxygen-documented; all public
  functions have proper `@param` / `@return` / `@details` / `@export`
  blocks.
- **Top-level orchestrator** `bisrDE::run_pipeline(counts_path,
  samplesheet_path, outdir, runid, annotation, id_type, brs_ticket)` —
  end-to-end pipeline as a callable function.
- **Quarto report template** at `bisrDE/inst/qmd/report.qmd` (with a
  `_sections/_comparison.qmd` child for the per-comparison loop). Replaces
  the parent's glue-built Rmd; uses native `params:`, `knit_child`, and
  Quarto callouts. Rendered via `bisrDE::generate_report()`.
- **Executive Summary** at the top of the report: programmatic 1-paragraph
  brief per comparison (DEG counts + top 3 up + top 3 down + strongest
  enriched gene set across the 4 backends).
- **Four-part interpretation callouts** for every plot section: "How to
  read this", "What to look for", "Common pitfalls", "Next steps" — for
  PCA, sample-level QC, volcano, heatmaps, and enrichment dotplots.
- **End-of-report Glossary** with 22 plain-English term definitions
  (padj, log2FC, NES, vst, TMM, GSEA, GO/KEGG/Reactome/Hallmark, FDR,
  z-score, etc.).
- **Annotation expansion** — every DE CSV now carries Ensembl ID, Entrez
  ID, gene Symbol, and gene name regardless of input identifier type.
- **`--id-type {ensembl,entrez,symbol}` flag** for non-Ensembl input
  count matrices. Symbol/Entrez inputs map to Ensembl IDs internally so
  downstream code that hardcodes `ENSEMBL_ID` keeps working.
- **`--brs-ticket BRS-XXXX` flag** plumbed through the pipeline; rendered
  as a YAML `subtitle:` under the report title.
- **Three new enrichment backends**: KEGG pathways (`gseKEGG`), Reactome
  pathways (`ReactomePA::gsePathway`), MSigDB Hallmark gene sets
  (`clusterProfiler::GSEA` + `msigdbr` H collection). Each produces its
  own dotplot and CSV alongside the existing GO results.
- **Top-100 heatmap variant** alongside the all-significant heatmap, with
  gene-symbol row labels for visual interpretation.
- **Sample-level QC plots** (4): Spearman correlation heatmap on DE genes,
  vst-transformed Euclidean distance heatmap, library-size + detected-gene
  barplot, hierarchical clustering (Ward.D2) + log-CPM density.
- **CLI overhaul** — replaced bare `print()` / `cat()` with `cli::*` for
  structured terminal output (banners, progress bars, ETA, color status,
  structured errors). Per-run session log file under `<outdir>/logs/`
  captures every cli/print/cat output via tee'd `sink()`.
- **Inline plotly PCA** — both 2D and 3D interactive plots now embedded
  in the report HTML, alongside the existing standalone HTML files.
- **Dynamic Methods text** — tool versions read from
  `bisrDE/inst/extdata/upstream_versions.yml` (single source of truth);
  R package versions read live via `packageVersion()`. Genome assembly
  string ("GRCh38 human" / "GRCm39 mouse") populated from the
  `--annotation` flag.
- **Nextflow DSL2 wrapper** at `differential_expression/nf-module/`:
  `main.nf` + `modules/local/bisr_de.nf` + `nextflow.config` (3 profiles:
  `local`, `container`, `slurm`). Drops in downstream of
  `nf-core/rnaseq`'s merged-counts output. See `nf-module/README.md`.
- **`nf-test` suite** at `nf-module/tests/main.nf.test` validating the
  process end-to-end against the bundled `assets/example_*` mouse fixture.
- **`testthat` suite** at `bisrDE/tests/testthat/` with 4 test files
  covering the bug-fix surface area: contrast parsing, annotation,
  heatmap sample-subset, volcano color/threshold alignment.
  `devtools::test()` reports `52 PASS / 0 FAIL`.
- **`dge_analysis.def` install step** for `bisrDE` via
  `remotes::install_local('/opt/dge_project/bisrDE', upgrade = 'never')`,
  preserving renv-pinned dep versions inside the container.

### Fixed

- **Volcano plot color/threshold mismatch** — y-axis, dot color, threshold
  line, and top-genes label selection all now use the same significance
  field (`sig`, default `padj`). The red dashed line now aligns with the
  grey/colored boundary it was supposed to represent.
- **Heatmap sample-subset bug** — the heatmap matrix is now restricted to
  the comparison's experimental + control samples only, instead of
  showing every sample in the full samplesheet (which previously made
  multi-contrast reports unreadable).
- **Manuscript text species mismatch** — mouse runs no longer say
  "GRCh38"; the genome-assembly string is now driven by the
  `--annotation` flag.
- **Report references renumbered + alphabetized** by first-author
  surname; every `[n]` in the Methods text now resolves uniquely. Added
  edgeR (Robinson) and ReactomePA (Yu) citations that were missing despite
  the tools being used.
- **Ambiguous one-to-many ID mappings** — `annotate_results` now resolves
  these by keeping the first match and emits a `cli_alert_warning` so the
  user is aware. Pre-checks input rownames for duplicates and warns.

### Changed

- **`de.R` reduced 1880 → 118 lines (94% reduction)** — now a thin
  optparse driver that delegates to `bisrDE::run_pipeline` +
  `bisrDE::generate_report`. The wrapper has a smart `library(bisrDE)` →
  `devtools::load_all("bisrDE")` fallback for local-dev iteration without
  reinstalls.
- **Report format**: Rmd → Quarto. Native `params:`, callouts, child
  templates eliminate the prior glue-concatenation bug class. Self-contained
  HTML (`embed-resources: true`) bundles every plot, table, and widget
  into one file.
- **Pipeline structure**: pre-existing functions migrated into per-domain
  R files under `bisrDE/R/` (io, normalize, de, annotate, enrich,
  plots_volcano, plots_heatmap, plots_pca, plots_qc, pipeline, cli, utils,
  report). NSE modernized to `.data$col` syntax for R CMD check
  compatibility; `%>%` → `|>` (R 4.2+).


### Deferred to HPC deployment

The package, Nextflow module, and local-R workflow are all complete and
verified independently. The container rebuild is a deployment artefact,
not a release artefact.

- **`dge_analysis.sif` container rebuild + in-container smoke**
  : def file is rebuild-ready (Quarto CLI install +
  lockfile collapse), but the actual `bash
  build_container.sh` invocation requires Linux x86_64 + Apptainer.
  Will land on the VCU HPRC login node at deployment time, where build
  env matches deploy env.

### Acceptance highlights

- `devtools::check()` 0 errors / 0 warnings / 1 NOTE (Quarto runtime deps,
  intentional).
- `devtools::test()` 52 PASS / 0 FAIL.
- `lintr::lint_package()` 271 → 20 lints (-93%, 0 dead-code findings).
- `bash run_analysis.sh ...` end-to-end smoke produces a 7.76 MB Quarto
  HTML report with embedded plotly widgets, base64 images, DT tables,
  Quarto callouts, complete TOC, BRS subtitle.
- `nextflow run nf-module/main.nf --help` lists all 7 params.
- `nf-test test nf-module/tests/main.nf.test` PASSES end-to-end against
  the example fixture (~90 s).

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

- Container rebuild with Quarto CLI + retry of Glimma 4.6.
- Multi-factor designs and covariate support.
- Batch correction options (sva / ComBat).
- Single-cell RNA-seq adapter (separate package, not bisrDE).

---

## Version History

| Version | Date       | Description                                                                                  |
| ------- | ---------- | -------------------------------------------------------------------------------------------- |
| 1.5.0   | 2026-06-12 | Sample/group/contrast selection, DisplayName plot labels, charmbracelet (gum) interactive CLI |
| 1.4.0   | 2026-05-01 | Refactor to `bisrDE` R package + Quarto report + Nextflow DSL2 wrapper + 4-backend GSEA + QC |
| 1.3.0   | 2026-02-02 | Stable release with full DE pipeline                                                         |
