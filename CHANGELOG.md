# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.7.0] - 2026-09-03

### Added

- **Three run options, threaded through `de.R`, `run_interactive.sh`, the Go
  TUI and `bisrDE::run_pipeline()`** (launcher/TUI `--print-cmd` parity kept):
  - `--gsea-rank stat|log2fc` (env `BISR_GSEA_RANK`, default `stat`): pre-ranked
    GSEA now ranks every tested gene by the DESeq2 Wald statistic instead of the
    unshrunken log2 fold change, so noisy low-count genes with large fold
    changes no longer dominate the tails. `log2fc` restores the previous
    behaviour. **Enrichment tables change for every run that keeps the
    default.**
  - `--lfc-shrink apeglm|normal|none` (env `BISR_LFC_SHRINK`, default
    `apeglm`): `DESeq2::lfcShrink()` adds `log2FC_shrunken` / `lfcSE_shrunken`
    to every DE CSV and places volcano points at the shrunken value; DEG calls,
    counts and heatmaps still use the unshrunken `log2FoldChange`. Falls back to
    `normal` with a message when apeglm is not installed. apeglm 1.20.0 added to
    the renv library and `Suggests`.
  - `--independent-filtering` (env `BISR_INDEPENDENT_FILTERING=yes`, default
    off): exposes DESeq2's `independentFiltering`. Off remains the BISR default
    so existing results are unchanged; the setting is now stated in the Methods
    and recorded in the run JSON.
- **QC at a glance.** `qc_at_a_glance()` computes the library-size and
  detected-gene ranges, PC1/PC2 variance, a between/within group separation
  ratio (with mean silhouette when `cluster` is available) and outlier flags.
  Outliers are scored by the distance to the nearest replicate of the same
  group (median + 3 MAD and at least twice the typical distance), because a
  centroid rule cannot fire in a 3-vs-3 design: the outlier drags its own
  centroid. The report prints the verdict above the PCA; the RDS bundle and
  the run JSON store the numbers.
- **Provenance in every artifact.** The RDS bundle keeps the fitted
  `DESeqDataSet`, design formula, reference level and coefficient per
  comparison, and the run JSON summarises them along with the run options,
  the platform, and the versions of R, bisrDE, DESeq2, apeglm, edgeR,
  clusterProfiler, ReactomePA, msigdbr and ComplexHeatmap. The report ends with
  a version table and a collapsed `sessionInfo()`, and its Methods text reads
  versions from the run record rather than from the rendering machine.
- **Input gene names are kept.** `read_counts()` returns a `gene_meta`
  attribute (versioned ID + input `gene_name`); `annotate_results()` uses it to
  name genes the OrgDb cannot (`SYMBOL_SOURCE` = `orgdb` / `input`) and adds
  `ENSEMBL_ID_VERSIONED` so a row can be traced to the nf-core matrix. An input
  name that is only the accession (nf-core's fallback when the GTF has no name)
  is not used, so `SYMBOL` stays `NA` for those genes.
- Tests: 84 new (238 total) covering the threshold predicate, log-scale
  z-scores, the GSEA ranking vector, group-size derivation, the QC summary
  including a 3-vs-3 outlier, the DESeq2 options, shrinkage and its failure
  path, coefficient naming, the QC matrix and its fallback, the volcano's
  shrunken axis, and the annotation fallback for Ensembl and symbol input. A
  bats case checks launcher/TUI parity with the three new flags set.

### Changed

- **Sample exploration runs on a log-scale, all-gene matrix.** PCA (static, 2D
  and 3D) and the Spearman correlation heatmap use the blind VST of the
  pre-filtered counts (`qc_matrix()`, falling back to `log2(TMM-CPM + 1)` on
  tiny inputs) instead of linear TMM-CPM; the correlation heatmap uses all
  genes, not the union of DE genes, which separated the groups by
  construction. The TMM CSV export is unchanged. Figures change; no DE
  statistic does.
- **DE heatmaps z-score `log2(TMM-CPM + 1)`** rather than jittered linear CPM;
  a zero-variance gene gets z = 0 instead of `NaN`.
- **Pre-filter group size follows the design**: `>= 10 reads in >= k samples`
  where `k` is the smallest group size (floor 2), previously a fixed 3.
- **GSEA p-values are exact** (`eps = 0` in every clusterProfiler / ReactomePA
  call), removing the 1e-10 floor that tied several Hallmark sets at
  p.adjust 1e-09. Dotplots, "top" tables and the Executive Summary now include
  only sets with `p.adjust < --padj`, ties on p.adjust break by |NES|, the
  summary table reports the number of significant sets per backend instead of
  Yes/No, and every enrichment CSV carries `significant` and `ranking_metric`
  columns.
- **Methods text is generated from the run**, not from a template: the
  smallest-group-size filter, count rounding, independent-filtering setting,
  shrinkage method, annotation collapse rule, GSEA metric and parameters, the
  QC transform and the platform are all read from the RDS bundle. The
  platform is recorded at analysis time by `.run_platform()`: the
  `BISR_PLATFORM_NAME` environment variable when set (export it in your job
  script, e.g. `"VCU's High Performance Research Computing cluster"`), else a
  Slurm description when `SLURM_JOB_ID` is set, else `Sys.info()` + R. The
  Pipeline section now says DESeq2 uses its own size factors and TMM is for
  plots and export.
- `filter_significant()` uses `padj < cutoff` like every other path (was
  `<=`), so a gene exactly at the cutoff is treated the same everywhere.
- `perform_deseq2_analysis()` returns the fitted `dds` and the options used as
  attributes; `run_analysis()` returns `dds` and `de_options`; `create_dotplot()`
  returns `NULL` when nothing is significant.

### Fixed

- **"Top 20 by padj" tables were mis-ordered.** They sorted the *formatted*
  padj string lexicographically, so `1.00e-02` came before `3.20e-35` and the
  strongest hits (e.g. CCND1, CDKN2A on the GSE72536 demo) were missing from
  the tables. They now sort numerically, then format.
- **Blank versions and `bisrDE vdev` in the Methods.** Two causes. (1)
  `upstream_versions.yml` carried a non-ASCII em dash in a comment; under the
  C locale of `run_analysis.sh` / Quarto subprocesses `yaml::read_yaml()`
  failed with "invalid input" and returned `NULL` silently, so every upstream
  version rendered as a blank (`STAR v [2]`). The file is ASCII-only now and
  is read as UTF-8 explicitly. (2) The render runs in its own R process where,
  in dev mode, `packageVersion("bisrDE")` and `system.file()` find nothing, so
  the package version rendered as `dev`. `generate_report()` now resolves the
  bisrDE version and the YAML path and passes them as params. A version that
  is missing renders as "version not recorded", never as a blank.
- **Hard-coded "performed on VCU's High Performance Research Computing
  cluster"** was asserted for every run, including laptop runs.
- Executive Summary no longer names the same comparison as both highest and
  lowest DEG count on single-contrast runs; the output tree is built from the
  directories that exist for the run and now lists the KEGG / Reactome /
  Hallmark / QC outputs.
- The one-to-many annotation warning said "ortholog"; it now says what it is.
- Per-comparison child sections are written with `writeLines(useBytes = TRUE)`
  and the two `≈` in the volcano callouts are HTML entities, so a report
  rendered under a C locale (Slurm, containers) no longer shows `<U+2248>`.

- **MSigDB Hallmark enrichment now works on offline HPC nodes.** `msigdbr` ≥ 24
  downloads its gene-set archive from Zenodo on first use, so on a compute node
  with no outbound internet every run failed with
  `Timeout was reached [zenodo.org]` and silently produced no Hallmark results.
  The archive is cached (`tools::R_user_dir("msigdbr", "cache")`) and the
  download is skipped entirely once present, so warming the cache once on a
  login node fixes it permanently (verified with outbound network fully
  blocked: Hallmark returns all 50 gene sets). Adds `just msigdb-cache` to warm
  it, README instructions (including `R_USER_CACHE_DIR` when `$HOME` is not
  shared with compute nodes), and, when the failure looks network related, a
  hint naming the cache path and the command to run instead of an opaque
  timeout.

## [1.6.4] - 2026-08-13

### Fixed

- **HPC install instructions no longer fail on the renv autoloader.** The
  documented `Rscript -e 'remotes::install_local("bisrDE")'` step failed on a
  fresh clone with `there is no package called 'remotes'`, even inside a
  correctly built conda environment: the shipped `.Rprofile` autoloads renv,
  which repoints `.libPaths()` at an empty project library and hides everything
  conda installed. `run_analysis.sh` already guarded pipeline *runs*, but not
  direct `Rscript` calls — which is exactly what the install step is. The README
  now writes the setting into the environment itself
  (`$CONDA_PREFIX/etc/conda/activate.d/renv_off.sh`, works for conda and
  micromamba), so every activation is safe, and documents the one-off
  `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` escape hatch alongside the exact error
  it resolves.

## [1.6.3] - 2026-08-13

HPC deployment: a conda environment, a validation gate, and a documented path
onto the cluster.

### Added

- **`environment.yml`** — a solve-verified conda environment for a centrally
  managed bioconda install (344 packages, verified with a real
  `micromamba --dry-run` against `linux-64`). It pins **Bioconductor 3.18 /
  R 4.3**, *not* `renv.lock`'s 3.16: Bioc 3.16 is **not installable** from
  bioconda, for two independent and unfixable reasons —
  `enrichplot 1.18.0` → `r-ggraph` → `r-ggforce` → `r-tweenr` has no
  R-4.2-compatible path, and `enrichplot` pulls `bioconductor-hdo.db`, which
  requires `annotationdbi >= 1.62`, newer than 3.16's 1.60. 3.18 is the nearest
  solvable release. The file documents the admin-facing caveats (named env not
  `base`, read-only, channel priority without `defaults`, `R_LIBS_USER`
  isolation).
- **`compare_runs.R`** — validation gate for an environment change. Diffs two
  run directories: gene coverage, log2FC Spearman/Pearson plus median and max
  drift, padj rank correlation, DEG-set Jaccard with the specific genes
  gained/lost, and enrichment term agreement across all four backends. Exits
  `0` PASS / `1` WARN / `2` FAIL so it can gate a deployment. Base R only, so
  it runs in either environment without becoming part of what it validates.
- **`submit_slurm.sh`** — sbatch template (4 cpu / 32 G / 4 h, overridable via
  `--export`) that handles the conda shell hook `conda activate` needs in a
  non-interactive job.
- README section documenting the full HPC path: clone → micromamba bootstrap
  (no admin rights) → environment → install `bisrDE` → smoke test → sbatch →
  validate.
- `just conda-env`, `just conda-solve <platform>`, `just compare <a> <b>`.

### Fixed

- **`run_analysis.sh` is now conda-aware.** The repo ships a `.Rprofile` +
  `renv/activate.R`, so on a fresh clone renv autoloaded and **shadowed an
  activated conda environment's R library** — producing "package not found"
  even when the environment was built correctly. R coming from `$CONDA_PREFIX`
  now disables the renv autoloader automatically, warns when `$CONDA_PREFIX` is
  set but `Rscript` resolves elsewhere, and skips the renv-bootstrap branch.
- `jsonlite` was used by `events.R` and `cli.R` (the `.report.json` summary and
  the progress event stream) but never declared in `DESCRIPTION`. It only
  worked because jsonlite arrived transitively; in a minimal environment both
  features would have broken.

## [1.6.2] - 2026-08-12

Launcher usability: know what the pipeline is doing, and stop typing paths.

### Added

- **Both launchers now find your input files for you.** They ask for a *project
  directory* first (default: cwd / the pipeline dir, or `BISR_BASE_DIR`), scan
  it, and offer the
  counts matrix and samplesheet as pick-lists instead of making you type full
  paths — best guess first, ranked by filename. `→ enter a path manually…` is
  always offered, and a directory with no candidates falls back to a plain text
  prompt. Crucially the scan **skips pipeline output** (`de_data`, `gsea_data`,
  `figures`, `logs`, `renv`, `work`, dotfiles), so a previous run's
  `DESeq2_*.csv` files can't swamp the samplesheet list. The bash launcher
  (`gum choose`) and the Go TUI use the same heuristics and agree on the same
  files — asserted by the bats suite. Interactive-only, so `--print-cmd` parity
  between the two front-ends is unaffected.

- **Phase status line under the progress bar.** The TUI now shows the current
  pipeline stage — `loading data`, `parsing contrasts`,
  `normalizing counts (TMM)`, `setting up DESeq2`,
  `running differential expression`, `generating QC plots`,
  `generating PCA plots`, `saving results`, `rendering report` — updating live
  as the run progresses. The spinner line names the *unit of work* (the current
  comparison) and the phase line names the *stage*, so nothing is duplicated.
  `generate_report()` emits its own phase, so the trailing Quarto render is
  visible with the bar at 100% instead of looking stalled.

## [1.6.1] - 2026-08-12

Live terminal UI, developer tooling, and real-data fixes.

### Added

- **Bounded live output box.** Pipeline output no longer scrolls the terminal
  forever: the most recent lines (default 8, `BISR_TUI_LINES`, clamped 3–40)
  are kept in a ring buffer and redrawn **in place** inside a bordered box that
  sits above the progress bar, with older lines falling off the top. The box
  keeps a stable height so the layout never jumps, and long lines are truncated
  to the terminal width. The terminal is never scrolled: lines that leave the
  top of the box are simply dropped from the view, and the complete transcript
  is written to `<outdir>/logs/<ts>_session.log` (plus `_session.html` when
  `aha` is installed).
- **Typewriter output streaming in the Go TUI.** Relayed pipeline lines are
  revealed a few characters at a time (~12 ms/char) instead of being dumped as
  chunks, so output flows organically. Tunable with `BISR_TUI_TYPE_DELAY_MS`
  (`0` disables it). The reveal rate **scales with backlog** and flushes
  outright past 40 queued lines, so a burst (Quarto alone emits hundreds of
  lines) can never leave the UI narrating output after the run has finished.
- **Functional severity colours.** Each line is classified by its `cli` glyph
  — header / success / info / warning / error — and given a deliberate style,
  so severity reads at a glance. Incoming ANSI is stripped first, which is what
  makes character-by-character streaming safe: revealing a raw escape sequence
  one byte at a time would otherwise spray garbage.
- **`bin/bisrde` wrapper + `just install-cli`** — `just` only finds its
  justfile by searching *upward*, so it works anywhere inside the repo but not
  outside it. The wrapper resolves the repo from its own location, so
  `bisrde tui` works from any directory.

- **Live event-driven Go TUI.** The pipeline now emits a structured NDJSON
  progress stream (`start` / `phase` / `tick` / `done`) when
  `BISR_EVENTS_FILE` is set, and the Go TUI tails it to drive a bubbletea UI:
  a `bubbles/progress` bar, spinner, current comparison and current step stay
  pinned beneath the live output box (see above). The UI is rendered from
  *state* rather than scraped from
  terminal output — the same model React/Ink gives Claude Code. When the event
  stream is active the R side suppresses its own `cli` bar, so exactly one
  component renders progress. Entirely opt-in: with the variable unset the
  bash launcher and `Rscript de.R` behave exactly as before, and any failure to
  start the UI falls back to plain output relaying.

- **`justfile`** — one discoverable surface (`just --list`) over the repo's
  entry points: `run`, `tui`, `demo`, `test`/`test-r`/`test-go`/`test-sh`,
  `check`, `doc`, `lint`, `container`, `bump`, plus gated recipes for the
  optional tooling (`inspect` visidata, `log-html` aha, `disk` dust/ncdu,
  `shot` pageres, `optimize` optimizt, `vhs`). A missing tool prints an install
  hint instead of failing obscurely.
- **bats-core parity suite** (`differential_expression/tests/bats/parity.bats`,
  `just test-sh`) — asserts `run_interactive.sh --print-cmd` and
  `bisrde-tui --print-cmd` emit **byte-identical** argument vectors across BRS
  ticket, threshold, volcano-label and contrast-selection combinations, that
  `--print-cmd` never executes the pipeline, and that shellcheck stays clean.
  This is the first automated guard on the two-launchers-drifting risk.
- **Colour-preserved HTML session log** — `stop_session_log()` renders
  `<outdir>/logs/<ts>_session.html` via `aha` from the raw log *before* ANSI is
  stripped, so you get both a highlighted, linkable page and a clean plain-text
  `.log`. No-op when `aha` is absent; failures never break a finished run.
- **`docs/demo.tape`** — a charmbracelet **VHS** script so the terminal demo
  re-records deterministically each release (`just vhs`) instead of being
  hand-performed.
- **Graceful interrupt** — Ctrl-C during a run now closes the session-log sink,
  finalizes the log (ANSI/`\r` stripped), reports where partial results were
  kept, and exits **130** instead of looking like a clean finish.

### Fixed

- **The comparisons progress bar is now genuinely live.** It advanced only once
  per comparison, so it froze for the minutes each comparison takes — cli never
  redraws on its own, a bar only repaints inside `cli_progress_update()`. Worse,
  with `total = n_comparisons` cli's default `auto_terminate` swallowed the
  final update, so a 3-contrast run only ever rendered 33% and 67% before
  "jumping" to done. `run_analysis()` now reports its 10 phases (DESeq2,
  annotate, save, volcano, 2 heatmaps, 4 enrichment backends) through a `tick`
  callback and the bar is sized `n_comparisons * 10`, with an initial forced 0%
  frame and `auto_terminate = FALSE`. Measured on the example dataset: in-place
  redraws 4 → 407, rendered frames 2 → 44, and the bar now names the running
  step. The session log stays ANSI/`\r`-free.

- **Samplesheet header variants no longer cause a misleading "No contrasts to
  analyse" failure.** A header like `Sample_ID,Group_ID` parsed without
  complaint but left `samplesheet$GroupID` `NULL`, so every contrast silently
  resolved no experimental/control group and the run died several steps later
  with an unrelated-sounding error. `read_samplesheet()` now normalizes common
  variants (`Sample_ID`, `sample id`, `Group_ID`, `condition`, …) to the
  canonical `SampleID` / `GroupID` / `DisplayName` / `Exclude` names — contrast
  columns (`*_vs_*`) are never renamed and no duplicate column can be created —
  and **aborts immediately with a clear message** if `SampleID`/`GroupID` are
  genuinely absent.
- The "No contrasts to analyse" abort is now self-diagnosing: it distinguishes
  "the samplesheet has no contrast columns" from "contrast columns found, but
  none resolved both an experimental and a control group", and says how to fix
  each.

## [1.5.4] - 2026-07-14

Addresses the first-run feedback in issue #12.

### Added

- **Configurable significance thresholds** — `--fold-change N` (linear, e.g. `2`
  for 2-fold; default `1.5`) and `--padj P` (default `0.05`) set the DEG cutoffs
  and are applied consistently to the volcano plots, heatmaps, DEG counts,
  `.report.json`, and the report's Methods text + callouts. Prompts added to
  both launchers (env `BISR_FOLD_CHANGE` / `BISR_PADJ`). (#12)

### Fixed

- `read_counts()` drops Ensembl `_PAR_Y` pseudo-autosomal genes, which otherwise
  collapse onto their chrX copy once the version suffix is stripped and errored
  with `duplicate 'row.names' are not allowed`. (#12)

### Changed

- **Flexible counts-file handling** — `read_counts()` now drops known metadata
  columns *by name* (`gene_name`, `transcript_id(s)`, `Chr`/`Start`/`End`/
  `Strand`/`Length`, …) so salmon / RSEM / featureCounts-style `nf-core/rnaseq`
  outputs all load without reformatting; the numeric `Start`/`End`/`Length`
  annotation columns are no longer mistaken for sample counts. The README
  documents the accepted layouts. (#12)

## [1.5.3] - 2026-07-07

### Added

- **Per-group display-name prompt** in both launchers (`#2b`). When the
  samplesheet has no usable `DisplayName` column, `run_interactive.sh` and the
  Go TUI offer **auto-derive** (default, `<GroupID> <n>`) or **enter a label per
  group** (per-sample doesn't scale) — the latter writes a working-copy
  samplesheet with `DisplayName = "<label> <n>"` and recommends adding the
  column for full control. Interactive-only, so `--print-cmd` parity is
  preserved and byte-identical across both front-ends.
- **Configurable volcano labels** — `--volcano-labels N` (default 10) caps how
  many genes are named on each volcano plot (top N by significance), with a
  prompt in `run_interactive.sh` and the Go TUI (env `BISR_VOLCANO_LABELS`).
- **Normalized counts in the DE spreadsheet** — each
  `DESeq2_<comparison>.csv` now carries per-sample normalized counts for all
  samples next to the log2FC / padj: `TMM_<SampleID>` (edgeR TMM) and
  `DESeq2norm_<SampleID>` (DESeq2 median-of-ratios).

### Fixed

- Volcano plots labelled *every* significant gene (and, for non-ensembl inputs,
  every unmapped gene via `NA %in% c(.., NA)`), burying the plot under hundreds
  of labels. They now label only the top `--volcano-labels` most-significant
  genes (default 10).

- `annotate_results()` now always emits a canonical `SYMBOL` (and `ENTREZID`)
  column regardless of `--id-type`. With `--id-type symbol` it previously
  produced only `SYMBOL_ID`, so the volcano / heatmap code (which uses
  `.data$SYMBOL`) errored and *every* comparison was caught and dropped as "no
  results" — leaving the report with no DE tables, GSEA, volcano plots, or
  dotplots. Symbol-input runs now complete end-to-end.
- The sample correlation heatmap and per-comparison DE heatmaps no longer
  silently skip for `--id-type symbol` / `entrez`: the significant-gene and
  heatmap subsets now match the count matrix by whichever DE ID column overlaps
  it (new `.match_id_column()` helper) instead of a hardcoded `ENSEMBL_ID`.
- `read_counts()` no longer collapses gene symbols that contain dots (mouse
  `H2-M10.1`, `Tex19.1`, `Rn4.5s`, … are *distinct* genes) into duplicate row
  names — the Ensembl version-suffix strip now applies only to `ENS…`
  accessions. Fixes a hard `duplicate 'row.names' are not allowed` error on
  `salmon.merged.gene_counts.tsv` files whose `gene_id` column holds symbols.
- Go TUI startup banner showed a stale `v1.5.0` (missed by the v1.5.2 version
  bump); now reads the current version.

## [1.5.2] - 2026-07-01

Report & UX fixes surfaced after v1.5.1.

### Fixed

- **Report figures were missing from the HTML.** `run_pipeline()` now
  normalizes `outdir` to an absolute path before building the figure
  directories, so the paths stored in the analysis RDS resolve when the Quarto
  report renders from a temp working dir. Previously every figure silently fell
  back to "not available".
- **Session log came out empty.** `cli` routes its output to `stderr` whenever
  a sink is active (NEWS #153), so the output-only session sink captured
  nothing. `start_session_log()` now pins cli to stdout via
  `start_app(output = "stdout", .auto_close = FALSE)` so narration is both shown
  live and written to the log; `stop_session_log()` strips ANSI for a clean
  plain-text file. The log now spans the report render too, not just analysis.

### Added

- **Auto-derived sample labels.** With no `DisplayName` column, plots now label
  samples `"<GroupID> <n>"` (a warning recommends adding the column) instead of
  showing raw accession SampleIDs. An explicit `DisplayName` still overrides per
  sample.
- **Colorblind-safe palette across all plots** — shared Okabe-Ito (categorical)
  and diverging / `Mako` (continuous) helpers applied to the QC plots, volcano,
  PCA (static + interactive), DE heatmap, and enrichment dotplots.
- **Structured run summary** at `<outdir>/logs/<ts>_report.json` — run metadata,
  timings, gene/sample counts, per-contrast DE (`n_sig`) + enrichment outcomes,
  and output paths.
- **Live Go TUI progress** — the TUI now runs the pipeline under a
  pseudo-terminal (`creack/pty`) so R's `cli` renders live progress bars instead
  of one-shot milestone prints (67% → 100%).
- `run_interactive.sh` screenshots the run-config card via charmbracelet
  `freeze` when installed (optional, gated, `--print-cmd` parity preserved).

### Changed

- QC plots widened; the vst-distance heatmap title no longer clips on the left
  edge and the hclust / density panel titles no longer overlap.

## [1.5.1] - 2026-06-12

### Added

- **Go TUI launcher** (charmbracelet Phase B) at `differential_expression/tui/`
  — a [bubbletea](https://github.com/charmbracelet/bubbletea)/[huh](https://github.com/charmbracelet/huh)/[lipgloss](https://github.com/charmbracelet/lipgloss)
  `bisrde-tui` binary with a styled form and interactive group / sample /
  contrast selection. Its assembled command is byte-identical to
  `run_interactive.sh --print-cmd` (parity is a maintained contract).
  `run_interactive.sh` now opens with a **bash-vs-Go-TUI chooser** that builds
  the TUI on demand and falls back to bash. Binary git-ignored; source committed.

### Changed

- `parse_contrasts()` warns when a contrast column has non-`0`/`1` values
  (previously a silent "missing exp or ctrl group").
- `run_pipeline()` warns on duplicate `DisplayName` labels.
- `run_interactive.sh` honors `BISR_INCLUDE_CONTRASTS` / `BISR_EXCLUDE_CONTRASTS`
  directly in non-interactive mode (parity with the Go TUI).

### Fixed

- Generated the roxygen man pages + the `filter_samplesheet` NAMESPACE export
  that were missing from v1.5.0 (committed without `devtools::document()`).

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

- Container rebuild + in-container smoke test on HPC (x86_64 / Apptainer).
- Multi-factor designs and covariate support.
- Batch correction options (sva / ComBat).
- Single-cell RNA-seq adapter (separate package, not bisrDE).

---

## Version History

| Version | Date       | Description                                                                                  |
| ------- | ---------- | -------------------------------------------------------------------------------------------- |
| 1.7.0   | 2026-09-03 | GSEA ranked by the Wald statistic with exact p-values, apeglm shrinkage, `--independent-filtering`; sample QC on blind VST (all genes); run provenance (dds, options, versions, QC summary) in RDS/JSON; report fixes (top-20 order, blank versions, Methods generated from the run, sessionInfo) |
| 1.6.4   | 2026-08-13 | Fix the HPC install step that renv's autoloader breaks (bake `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` into the conda env) |
| 1.6.3   | 2026-08-13 | HPC deployment: solve-verified conda `environment.yml` (Bioc 3.18), `compare_runs.R` validation gate, Slurm template, conda-aware launcher |
| 1.6.2   | 2026-08-12 | Phase status line under the progress bar; both launchers discover counts/samplesheet files from a project directory instead of requiring typed paths |
| 1.6.1   | 2026-08-12 | Live event-driven Go TUI (NDJSON progress → bubbletea box + bar, typewriter streaming, severity colours); justfile / bats parity suite / shellcheck / aha / VHS; graceful interrupt; samplesheet header variants; genuinely live progress |
| 1.5.4   | 2026-07-14 | Configurable DE thresholds (`--fold-change`/`--padj`, dynamic throughout the report); flexible counts-file layouts (salmon/RSEM/featureCounts); drop Ensembl `_PAR_Y` (#12) |
| 1.5.3   | 2026-07-07 | `--id-type symbol`/`entrez` fixes (SYMBOL col, heatmaps, dotted-symbol counts); volcano top-N labels (`--volcano-labels`); per-group DisplayName prompt; normalized counts (TMM + DESeq2) in DE sheets |
| 1.5.2   | 2026-07-01 | Report figures embed (abs outdir); auto-derived + colorblind-safe plots; live Go TUI (PTY); non-empty session log + `.report.json` |
| 1.5.1   | 2026-06-12 | Go TUI launcher (charmbracelet Phase B); parse_contrasts / DisplayName warnings; roxygen + NAMESPACE fixes |
| 1.5.0   | 2026-06-12 | Sample/group/contrast selection, DisplayName plot labels, charmbracelet (gum) interactive CLI |
| 1.4.0   | 2026-05-01 | Refactor to `bisrDE` R package + Quarto report + Nextflow DSL2 wrapper + 4-backend GSEA + QC |
| 1.3.0   | 2026-02-02 | Stable release with full DE pipeline                                                         |
