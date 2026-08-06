# Top-level pipeline orchestrator: directory setup, per-comparison driver,
# and the end-to-end `run_pipeline()` that wraps the parent's main block.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.5 spec.
#
# Phase 5.5 design changes vs. parent:
# - `run_analysis()` takes `annotation`, `annotation_db`, `id_type` as
#   explicit args (parent inherited them as globals). The package version
#   threads them through to `annotate_results()` and the enrich fns so
#   the call chain has no hidden state.
# - `run_pipeline()` is NEW — wraps the parent's main block as a callable
#   function. Stops at the RDS save; report rendering is Phase 5.6 scope.

#' Set up the analysis output directory tree
#'
#' @description Create the standard set of subdirectories (data, figures,
#'   logs) under `base_dir`. Idempotent — existing directories are left
#'   alone.
#'
#' @param base_dir Base directory path. Sub-paths are created with
#'   `recursive = TRUE`.
#' @return A named list of created directory paths: `data`, `figures`,
#'   `de_data`, `gsea_data`, `kegg_data`, `reactome_data`,
#'   `hallmark_data`, `volcano`, `heatmap`, `gsea`, `kegg`, `reactome`,
#'   `hallmark`, `pca`, `qc`, `logs`.
#' @export
setup_directories <- function(base_dir) {
  dirs <- list(
    data           = file.path(base_dir, "data"),
    figures        = file.path(base_dir, "figures"),
    de_data        = file.path(base_dir, "data/de_data"),
    gsea_data      = file.path(base_dir, "data/gsea_data"),
    kegg_data      = file.path(base_dir, "data/kegg_data"),
    reactome_data  = file.path(base_dir, "data/reactome_data"),
    hallmark_data  = file.path(base_dir, "data/hallmark_data"),
    volcano        = file.path(base_dir, "figures/volcano"),
    heatmap        = file.path(base_dir, "figures/heatmap"),
    gsea           = file.path(base_dir, "figures/gsea"),
    kegg           = file.path(base_dir, "figures/kegg"),
    reactome       = file.path(base_dir, "figures/reactome"),
    hallmark       = file.path(base_dir, "figures/hallmark"),
    pca            = file.path(base_dir, "figures/pca"),
    qc             = file.path(base_dir, "figures/qc"),
    logs           = file.path(base_dir, "logs")
  )

  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}


#' Number of progress sub-steps [run_analysis()] reports per comparison
#'
#' @description Keep in sync with the `.tick()` calls in [run_analysis()]:
#'   DESeq2, annotate, save DE table, volcano, 2 heatmaps, and 4 enrichment
#'   backends. Used by [run_pipeline()] to size the progress bar so it advances
#'   during a comparison rather than once per comparison.
#' @keywords internal
.STEPS_PER_COMPARISON <- 10L


#' Run the full per-comparison DE + enrichment workflow
#'
#' @description Per-comparison driver: DESeq2 contrast, annotation, volcano
#'   plot, two heatmaps (all-significant + top-100), and four enrichment
#'   analyses (GO, KEGG, Reactome, MSigDB Hallmark) each with a dotplot.
#'   Catches all errors at the comparison boundary and returns `NULL` so
#'   one bad contrast doesn't take down the entire pipeline.
#'
#' @param comparison List with `name`, `exp`, `ctrl` (typically one entry
#'   from [parse_contrasts()] output).
#' @param dds A `DESeqDataSet` (post pre-filtering, pre-`DESeq()`).
#' @param normalized_counts TMM-normalized count matrix (Genes x Samples).
#' @param sample_info Data frame with `sample` + `condition` columns.
#' @param out_dirs Output directory list from [setup_directories()].
#' @param annotation One of `"human"` or `"mouse"` — used for KEGG /
#'   Reactome / MSigDB organism dispatch.
#' @param annotation_db `OrgDb` object (e.g. `org.Hs.eg.db::org.Hs.eg.db`)
#'   for `annotate_results` and `setReadable` Symbol mapping.
#' @param id_type Identifier type for the count rownames (one of
#'   `"ensembl"`, `"entrez"`, `"symbol"`). Default `"ensembl"`.
#' @return On success, a named list `list(deseq, gsea, kegg, reactome,
#'   hallmark)` where each enrichment slot is either a `setReadable`'d
#'   GSEA object or `NULL` if no enrichment was found. On error, returns
#'   `NULL` and emits `cli::cli_alert_danger`.
#'
#' @importFrom utils write.csv
#' @importFrom grDevices png dev.off
#' @export
run_analysis <- function(comparison, dds, normalized_counts, sample_info,
                         out_dirs, annotation, annotation_db,
                         id_type = "ensembl", volcano_labels = 10,
                         deseq_norm_counts = NULL,
                         padj = 0.05, lfc = 0.58,
                         tick = NULL) {
  # Advance the caller's progress bar one sub-step. No-op when called directly
  # (tick = NULL), so this function still works outside run_pipeline().
  .tick <- function(step) if (is.function(tick)) tick(step)

  tryCatch(
    {
      cli::cli_h2("Comparison: {.strong {comparison$name}}")
      cli::cli_inform(c(
        "*" = "Experimental: {.val {comparison$exp}}",
        "*" = "Control:      {.val {comparison$ctrl}}"
      ))

      .tick("DESeq2")
      deseq_results <- perform_deseq2_analysis(dds, comparison$exp, comparison$ctrl)
      if (is.null(deseq_results)) {
        cli::cli_alert_warning("DESeq2 analysis returned NULL; skipping comparison")
        return(NULL)
      }

      .tick("annotate")
      cli::cli_alert_info("Annotating results...")
      annotated_results <- annotate_results(
        deseq_results,
        id_type       = id_type,
        annotation_db = annotation_db
      )

      # DE spreadsheet = annotated results + per-sample normalized counts (TMM +
      # DESeq2 median-of-ratios, all samples), so expression sits alongside the
      # log2FC. Built on a copy so the plots below use the un-widened table.
      de_out  <- annotated_results
      key_col <- .match_id_column(de_out, rownames(normalized_counts))
      if (!is.null(key_col)) {
        ids      <- as.character(de_out[[key_col]])
        orig_ids <- stats::setNames(as.character(sample_info$sample),
                                    make.names(as.character(sample_info$sample)))
        de_out <- .append_counts(de_out, ids, normalized_counts, orig_ids, "TMM")
        if (!is.null(deseq_norm_counts)) {
          de_out <- .append_counts(de_out, ids, deseq_norm_counts, orig_ids, "DESeq2norm")
        }
      }

      output_file <- create_file_path(out_dirs$de_data, "DESeq2_", comparison$name)
      .tick("save DE table")
      cli::cli_alert_info("Saving DE results: {.path {output_file}}")
      utils::write.csv(de_out, output_file)

      .tick("volcano")
      cli::cli_alert_info("Generating volcano plot...")
      volcano_plot <- generate_volcano(
        annotated_results,
        comparison$exp,
        comparison$ctrl,
        p = padj, lfc = lfc,
        n_labels = volcano_labels
      )
      save_plot(
        volcano_plot,
        create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png")
      )

      .tick("heatmap (all)")
      cli::cli_alert_info("Generating heatmap (all significant DEGs)...")
      grDevices::png(
        create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap_all_sig.png"),
        width = 800, height = 1200, res = 150
      )
      generate_heatmap(annotated_results, normalized_counts, sample_info,
                       p = padj, lfc = lfc,
                       exp_name = comparison$exp, ctrl_name = comparison$ctrl)
      grDevices::dev.off()

      .tick("heatmap (top 100)")
      cli::cli_alert_info("Generating heatmap (top 100 DEGs by padj)...")
      grDevices::png(
        create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap_top100.png"),
        width = 800, height = 1400, res = 150
      )
      generate_heatmap(annotated_results, normalized_counts, sample_info,
                       p = padj, lfc = lfc,
                       exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                       top_n = 100)
      grDevices::dev.off()

      # ---- Enrichment: GO ----
      .tick("GO")
      cli::cli_alert_info("Enrichment: GO (gseGO)...")
      gse <- process_gsea(annotated_results, annotation_db = annotation_db)
      if (!is.null(gse)) {
        utils::write.csv(
          as.data.frame(gse),
          create_file_path(out_dirs$gsea_data, "GO_Analysis_", comparison$name)
        )
        cli::cli_alert_info("  Generating GO dotplot...")
        save_plot(
          create_dotplot(gse, create_comparison_name(comparison$exp,
                                                     comparison$ctrl, "GSEA-GO ")),
          create_file_path(out_dirs$gsea, "", comparison$name, "_GSEA.png")
        )
      } else {
        cli::cli_alert_warning("  Skipping GO dotplot - no enrichment results")
      }

      # ---- Enrichment: KEGG ----
      .tick("KEGG")
      cli::cli_alert_info("Enrichment: KEGG (gseKEGG)...")
      kegg_gse <- process_kegg_gsea(annotated_results,
                                    annotation = annotation,
                                    annotation_db = annotation_db)
      if (!is.null(kegg_gse)) {
        utils::write.csv(
          as.data.frame(kegg_gse),
          create_file_path(out_dirs$kegg_data, "KEGG_Analysis_", comparison$name)
        )
        cli::cli_alert_info("  Generating KEGG dotplot...")
        save_plot(
          create_dotplot(kegg_gse, create_comparison_name(comparison$exp,
                                                          comparison$ctrl, "KEGG ")),
          create_file_path(out_dirs$kegg, "", comparison$name, "_KEGG.png")
        )
      } else {
        cli::cli_alert_warning("  Skipping KEGG dotplot - no enrichment results")
      }

      # ---- Enrichment: Reactome ----
      .tick("Reactome")
      cli::cli_alert_info("Enrichment: Reactome (gsePathway)...")
      reactome_gse <- process_reactome_gsea(annotated_results,
                                            annotation = annotation,
                                            annotation_db = annotation_db)
      if (!is.null(reactome_gse)) {
        utils::write.csv(
          as.data.frame(reactome_gse),
          create_file_path(out_dirs$reactome_data, "Reactome_Analysis_",
                           comparison$name)
        )
        cli::cli_alert_info("  Generating Reactome dotplot...")
        save_plot(
          create_dotplot(reactome_gse,
                         create_comparison_name(comparison$exp, comparison$ctrl,
                                                "Reactome ")),
          create_file_path(out_dirs$reactome, "", comparison$name, "_Reactome.png")
        )
      } else {
        cli::cli_alert_warning("  Skipping Reactome dotplot - no enrichment results")
      }

      # ---- Enrichment: MSigDB Hallmark ----
      .tick("Hallmark")
      cli::cli_alert_info("Enrichment: MSigDB Hallmark...")
      hallmark_gse <- process_msigdb_hallmark(annotated_results,
                                              annotation = annotation,
                                              annotation_db = annotation_db)
      if (!is.null(hallmark_gse)) {
        utils::write.csv(
          as.data.frame(hallmark_gse),
          create_file_path(out_dirs$hallmark_data, "Hallmark_Analysis_",
                           comparison$name)
        )
        cli::cli_alert_info("  Generating Hallmark dotplot...")
        save_plot(
          create_dotplot(hallmark_gse,
                         create_comparison_name(comparison$exp, comparison$ctrl,
                                                "Hallmark ")),
          create_file_path(out_dirs$hallmark, "", comparison$name, "_Hallmark.png")
        )
      } else {
        cli::cli_alert_warning("  Skipping Hallmark dotplot - no enrichment results")
      }

      cli::cli_alert_success("Comparison {.strong {comparison$name}} complete")
      list(
        deseq    = annotated_results,
        gsea     = gse,
        kegg     = kegg_gse,
        reactome = reactome_gse,
        hallmark = hallmark_gse
      )
    },
    error = function(e) {
      cli::cli_alert_danger("run_analysis error: {e$message}")
      NULL
    }
  )
}


#' Top-level bulk RNA-seq DE pipeline
#'
#' @description End-to-end bulk RNA-seq differential expression pipeline:
#'   read counts + samplesheet, parse contrasts, TMM normalize, run DESeq2
#'   per comparison (with GO/KEGG/Reactome/MSigDB enrichment), produce
#'   sample-level QC plots and PCA, save an RDS bundle for downstream
#'   reporting.
#'
#' @param counts_path Path to the merged counts TSV (typically
#'   `rsem.merged.gene_counts.tsv` from nf-core/rnaseq).
#' @param samplesheet_path Path to the comma-delimited samplesheet.
#' @param outdir Base output directory. Subdirectories created via
#'   [setup_directories()].
#' @param runid Run identifier (used in cli output and downstream report
#'   metadata).
#' @param annotation One of `"human"` or `"mouse"`. Drives both organism
#'   dispatch in the enrich fns and the OrgDb selection.
#' @param id_type One of `"ensembl"`, `"entrez"`, `"symbol"`. Identifies
#'   the gene-ID type used as count rownames. Default `"ensembl"`.
#' @param brs_ticket Optional BRS ticket identifier (e.g. `"BRS-1234"`).
#'   Default `""`. Currently passed through into the returned artifact
#'   list for the report renderer.
#' @param exclude_samples Optional character vector of `SampleID`s to drop
#'   from the entire analysis (see [filter_samplesheet()]). Default `NULL`.
#' @param exclude_groups Optional character vector of `GroupID`s to drop.
#'   Default `NULL`.
#' @param include_contrasts Optional character vector. When non-NULL, only
#'   these contrast columns are processed (allowlist). Default `NULL`.
#' @param exclude_contrasts Optional character vector of contrast columns to
#'   skip (denylist). Default `NULL`.
#' @param volcano_labels Max genes to label on each volcano plot (the top N by
#'   significance). Default `10`. Passed through to [generate_volcano()].
#' @param tick Optional function of one argument (a step label) called before
#'   each of the [.STEPS_PER_COMPARISON] phases, so the caller can advance a
#'   progress bar mid-comparison. `NULL` (default) disables it.
#' @param padj Adjusted p-value (FDR) cutoff for calling DEGs. Default `0.05`.
#' @param fold_change Linear fold-change cutoff for calling DEGs (e.g. `2` for
#'   2-fold). Default `1.5`; converted to a log2 cutoff internally and applied
#'   to the volcano / heatmaps / DEG counts / report text.
#' @param session_log Optional session-log handle from [start_session_log()].
#'   When supplied, this call uses it (so a driver can make one log span both
#'   the analysis and the report render) and does NOT close it — the caller
#'   owns the lifecycle. When `NULL` (default) the pipeline opens and closes
#'   its own log, as before. Default `NULL`.
#' @return Invisibly, a named list of pipeline artifacts:
#'   `results`, `comparisons`, `out_dirs`, `pca_plot` (ggplot),
#'   `pca_plotly` (plotly 2D), `pca_plotly_3d` (plotly 3D),
#'   `annotation`, `runid`, `brs_ticket`, `sample_info`, `tmm`,
#'   `counts`, `dds`, `rds_path` (path to the saved analysis RDS).
#' @details Loads the appropriate OrgDb (`org.Hs.eg.db` or `org.Mm.eg.db`)
#'   via `requireNamespace` — both packages are listed in `Suggests:`,
#'   not `Imports:`, since each run only needs one. Stops with a clear
#'   message if the required OrgDb is not installed.
#'
#'   Stops at the RDS save. Report rendering is a separate concern (Phase
#'   5.6 will introduce `bisrDE::generate_report()`); for now, downstream
#'   callers handle report generation themselves.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom utils write.csv
#' @export
run_pipeline <- function(counts_path,
                         samplesheet_path,
                         outdir,
                         runid,
                         annotation,
                         id_type           = "ensembl",
                         brs_ticket        = "",
                         exclude_samples   = NULL,
                         exclude_groups    = NULL,
                         include_contrasts = NULL,
                         exclude_contrasts = NULL,
                         volcano_labels    = 10,
                         padj              = 0.05,
                         fold_change       = 1.5,
                         session_log       = NULL) {
  annotation <- match.arg(annotation, c("human", "mouse"))
  id_type    <- match.arg(id_type, c("ensembl", "entrez", "symbol"))
  run_started <- Sys.time()
  # Significance thresholds: `fold_change` is the linear cutoff (user-facing,
  # e.g. 2 for 2-fold); `lfc` is its log2 form used by the plots/filters.
  lfc <- log2(fold_change)

  # ---- 1. Output dirs + session log ----
  # outdir MUST be absolute: the figure paths derived from it are stored in
  # the analysis RDS and resolved later by generate_report(), which renders
  # from a fresh temp working dir. A relative outdir would make every figure
  # path fail file.exists() at render time, so the report would silently show
  # "not available" for plots that are in fact on disk.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir   <- normalizePath(outdir, mustWork = FALSE)
  out_dirs <- setup_directories(outdir)
  # The driver (de.R) may pass an already-open session log so it spans both
  # this call AND the report render; only manage our own when it doesn't.
  own_log <- is.null(session_log)
  if (own_log) {
    session_log <- start_session_log(out_dirs$logs)
    on.exit(stop_session_log(session_log), add = TRUE)
  }

  cli::cli_h1("Bulk RNA-Seq Differential Expression Pipeline")
  cli::cli_inform(c(
    "Run ID:        {.val {runid}}",
    "Counts file:   {.path {counts_path}}",
    "Samplesheet:   {.path {samplesheet_path}}",
    "Output dir:    {.path {outdir}}",
    "Annotation:    {.val {annotation}}",
    "ID type:       {.val {id_type}}",
    if (nzchar(brs_ticket)) "BRS ticket:    {.val {brs_ticket}}" else NULL
  ))
  cli::cli_alert_info("Session log: {.path {session_log$path}}")

  # ---- 2. Load OrgDb (Suggests-style) ----
  annotation_db <- if (annotation == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' is required for annotation = 'human'. ",
           "Install with: BiocManager::install('org.Hs.eg.db')")
    }
    get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
  } else {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package 'org.Mm.eg.db' is required for annotation = 'mouse'. ",
           "Install with: BiocManager::install('org.Mm.eg.db')")
    }
    get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
  }

  # ---- 3. Read data ----
  cli::cli_h1("Loading data")
  counts <- read_counts(counts_path)
  samplesheet <- read_samplesheet(samplesheet_path)
  cli::cli_alert_success(
    "Loaded {nrow(counts)} genes x {ncol(counts)} samples; {nrow(samplesheet)} samplesheet rows"
  )

  # ---- 3b. Apply sample / group exclusions (Exclude column + flags) ----
  if (length(exclude_samples) > 0 || length(exclude_groups) > 0 ||
      "Exclude" %in% colnames(samplesheet)) {
    cli::cli_h1("Applying sample exclusions")
    samplesheet <- filter_samplesheet(
      samplesheet,
      exclude_samples = exclude_samples,
      exclude_groups  = exclude_groups
    )
  }

  # ---- 4. Parse contrasts (with optional include/exclude filtering) ----
  cli::cli_h1("Parsing contrasts")
  comparisons <- parse_contrasts(
    samplesheet,
    include_contrasts = include_contrasts,
    exclude_contrasts = exclude_contrasts
  )
  if (length(comparisons) == 0) {
    # Say WHY, not just that it happened: the common causes are a samplesheet
    # with no contrast columns at all vs. contrast columns whose 1/0 coding
    # never yields both an experimental and a control group.
    n_contrast_cols <- length(setdiff(
      colnames(samplesheet)[3:ncol(samplesheet)],
      .samplesheet_meta_cols
    ))
    cli::cli_abort(c(
      "No contrasts to analyse after include/exclude filtering and group checks.",
      if (n_contrast_cols == 0) c(
        "x" = "The samplesheet has no contrast columns.",
        "i" = "Add one column per comparison, named {.val <experimental>_vs_<control>}, after {.val GroupID}."
      ) else c(
        "x" = "{n_contrast_cols} contrast column{?s} found, but none resolved BOTH an experimental and a control group.",
        "i" = "In each contrast column mark experimental samples {.val 1}, control samples {.val 0}, and leave non-participating samples blank.",
        "i" = "Check that {.val GroupID} is populated for the marked rows."
      )
    ))
  }

  # ---- 5. Coerce + align ----
  countsdf <- counts |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),
                                \(x) as.integer(round(x))))
  countsdf <- align_counts_to_samplesheet(countsdf, samplesheet)

  # ---- 6. TMM normalization ----
  cli::cli_h1("Normalization (edgeR TMM)")
  tmm <- run_tmm(countsdf)
  utils::write.csv(
    tmm,
    file.path(out_dirs$de_data,
              paste0("normalizedCounts_tmm", Sys.Date(), ".csv"))
  )
  cli::cli_alert_success("TMM-normalized counts saved")

  # ---- 7. DESeq2 dataset + pre-filter ----
  cli::cli_h1("DESeq2 setup")
  sample_info <- data.frame(
    sample    = samplesheet$SampleID,
    condition = samplesheet$GroupID,
    stringsAsFactors = FALSE
  )
  # Per-sample display labels for plots (matrix join key stays SampleID; only
  # the visible labels change). Default to friendly auto-derived
  # "<GroupID> <n>" labels so plots never show raw accession SampleIDs; an
  # explicit DisplayName column overrides per sample (blank/NA cells keep the
  # auto-derived label).
  disp <- .auto_display(sample_info$condition)
  if ("DisplayName" %in% colnames(samplesheet)) {
    d   <- trimws(as.character(samplesheet$DisplayName))
    has <- !is.na(d) & nzchar(d)
    disp[has] <- d[has]
    cli::cli_alert_info(
      "DisplayName column detected: using it for plot labels (blank cells auto-derived)."
    )
  } else {
    cli::cli_alert_warning(c(
      "No DisplayName column: plot labels auto-derived from GroupID + replicate (e.g. {.val {disp[1]}}).",
      "i" = "Add a DisplayName column to the samplesheet for custom labels."
    ))
  }
  sample_info$display <- disp
  dups <- unique(disp[duplicated(disp)])
  if (length(dups) > 0) {
    cli::cli_alert_warning(
      "Duplicate display label{?s} {.val {dups}}; plot labels will be ambiguous for those samples."
    )
  }

  # Defensive: counts columns and colData rows must line up 1:1. After
  # exclusions this is the place a desync would surface, so fail with a
  # clear message rather than DESeq2's opaque "ncol == nrow is not TRUE".
  if (ncol(countsdf) != nrow(sample_info)) {
    cli::cli_abort(c(
      "Count columns ({ncol(countsdf)}) do not match samplesheet rows ({nrow(sample_info)}).",
      "i" = "This usually means sample alignment / exclusion left the two out of sync."
    ))
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countsdf,
    colData   = sample_info,
    design    = ~condition
  )

  smallest_group_size <- 3
  keep <- rowSums(countsdf >= 10) >= smallest_group_size
  dds <- dds[keep, ]
  cli::cli_inform(c(
    "Pre-filter: kept {sum(keep)} / {length(keep)} genes (>= 10 reads in >= {smallest_group_size} samples)",
    "Condition levels: {.val {levels(dds$condition)}}"
  ))

  # DESeq2 median-of-ratios normalized counts for the DE spreadsheets. Size
  # factors are global, so compute once here and reuse across comparisons.
  deseq_norm_counts <- DESeq2::counts(DESeq2::estimateSizeFactors(dds),
                                      normalized = TRUE)

  # ---- 8. Per-comparison loop ----
  cli::cli_h1("Per-comparison differential expression")
  results <- vector("list", length(comparisons))

  # Progress granularity: cli NEVER redraws on its own — the bar only repaints
  # inside cli_progress_update(). Ticking once per comparison therefore froze
  # the bar for the minutes each comparison takes, and with total = n the
  # default `auto_terminate` swallowed the final update (3 comparisons only
  # ever rendered 33% and 67%, then "jumped" to done). Counting sub-steps
  # instead gives ~.STEPS_PER_COMPARISON repaints per comparison.
  pb_total <- length(comparisons) * .STEPS_PER_COMPARISON
  pb <- cli::cli_progress_bar(
    name   = "Comparisons",
    total  = pb_total,
    format = paste0(
      "{cli::pb_name} {cli::pb_extra$current} | {cli::pb_extra$step} | ",
      "[{cli::pb_bar}] {cli::pb_percent} | ETA {cli::pb_eta}"
    ),
    extra  = list(current = "", step = "starting"),
    clear  = FALSE,
    # Default TRUE would end (and destroy) the bar the moment the last tick
    # reaches `total`, so the end-of-comparison re-sync below would then fail
    # with "Cannot find progress bar". Terminate explicitly via
    # cli_progress_done() after the loop instead.
    auto_terminate = FALSE
  )
  # cli_progress_bar() paints nothing; without this the bar first appears
  # already part-way through. Force an initial 0% frame.
  cli::cli_progress_update(id = pb, set = 0, force = TRUE)

  for (i in seq_along(comparisons)) {
    .label <- sprintf("%d/%d %s", i, length(comparisons), comparisons[[i]]$name)
    tick <- function(step) {
      cli::cli_progress_update(
        id = pb, inc = 1,
        extra = list(current = .label, step = step)
      )
    }
    res <- run_analysis(
      comparison        = comparisons[[i]],
      dds               = dds,
      normalized_counts = tmm,
      sample_info       = sample_info,
      out_dirs          = out_dirs,
      annotation        = annotation,
      annotation_db     = annotation_db,
      id_type           = id_type,
      volcano_labels    = volcano_labels,
      deseq_norm_counts = deseq_norm_counts,
      padj              = padj,
      lfc               = lfc,
      tick              = tick
    )
    # Re-sync: a comparison that bailed early (e.g. DESeq2 returned NULL) will
    # have ticked fewer than .STEPS_PER_COMPARISON times, so pin the bar to the
    # exact position for i completed comparisons.
    cli::cli_progress_update(
      id = pb, set = i * .STEPS_PER_COMPARISON, force = TRUE,
      extra = list(current = .label, step = "done")
    )
    if (!is.null(res)) {
      results[[i]] <- res
      cli::cli_alert_success(paste0(
        "[{i}/{length(comparisons)}] {.strong {comparisons[[i]]$name}}: ",
        "{nrow(res$deseq)} genes; ",
        "GO={!is.null(res$gsea)} ",
        "KEGG={!is.null(res$kegg)} ",
        "Reactome={!is.null(res$reactome)} ",
        "Hallmark={!is.null(res$hallmark)}"
      ))
    } else {
      cli::cli_alert_warning(
        "[{i}/{length(comparisons)}] {.strong {comparisons[[i]]$name}}: no results"
      )
    }
  }
  cli::cli_progress_done()

  # ---- 9. Sample Exploration QC plots ----
  cli::cli_h1("Sample Exploration QC plots")

  # These are intersected against rownames(tmm) — the INPUT gene IDs — so match
  # on whichever DE ID column overlaps the count matrix (ensembl / symbol /
  # entrez per --id-type). A hardcoded ENSEMBL_ID never matches a symbol/entrez
  # count matrix, which silently skipped the correlation heatmap.
  sig_genes_union <- unique(unlist(lapply(results, function(r) {
    if (is.null(r) || is.null(r$deseq)) return(character(0))
    d <- r$deseq
    ok <- !is.na(d$padj) & d$padj < padj & abs(d$log2FoldChange) >= lfc
    if (!any(ok)) return(character(0))
    key_col <- .match_id_column(d, rownames(tmm))
    if (!is.null(key_col)) as.character(d[[key_col]][ok]) else rownames(d)[ok]
  })))
  if (length(sig_genes_union) < 2) {
    cli::cli_alert_info(
      "No significant genes union - using all genes in TMM matrix for correlation heatmap"
    )
    sig_genes_union <- rownames(tmm)
  }
  cli::cli_alert_info(
    "Correlation heatmap will use {length(sig_genes_union)} gene{?s}"
  )

  qc_correlation_heatmap(
    normalized_counts = tmm,
    sample_info       = sample_info,
    sig_genes         = sig_genes_union,
    fig_path          = file.path(out_dirs$qc, "qc_correlation_heatmap.png")
  )
  cli::cli_alert_success("QC: correlation heatmap saved")

  qc_vst_dist_heatmap(
    dds         = dds,
    sample_info = sample_info,
    fig_path    = file.path(out_dirs$qc, "qc_vst_dist_heatmap.png")
  )
  cli::cli_alert_success("QC: vst distance heatmap saved")

  qc_libsize_detected_barplot(
    counts      = countsdf,
    sample_info = sample_info,
    fig_path    = file.path(out_dirs$qc, "qc_libsize_detected.png")
  )
  cli::cli_alert_success("QC: library size + detected genes barplot saved")

  qc_hclust_density(
    tmm         = tmm,
    sample_info = sample_info,
    fig_path    = file.path(out_dirs$qc, "qc_hclust_density.png")
  )
  cli::cli_alert_success("QC: hclust + density saved")

  # ---- 10. PCA ----
  cli::cli_h1("PCA")
  pca_plot <- pca_static(tmm, sample_info)
  save_plot(pca_plot, file.path(out_dirs$pca, "PCA_plot.png"))
  cli::cli_alert_success("Static 2D PCA saved")

  pca_plotly_2d_obj <- pca_plotly(tmm, sample_info)
  pca_plotly_2d_path <- file.path(out_dirs$pca, "allsamples_PCA_plot.html")
  export_plotly_to_html(pca_plotly_2d_obj, pca_plotly_2d_path)
  cli::cli_alert_success(
    "Interactive 2D PCA saved: {.path {pca_plotly_2d_path}}"
  )

  pca_plotly_3d_obj <- pca_plotly_3d(tmm, sample_info)
  pca_plotly_3d_path <- file.path(out_dirs$pca, "allsamples_PCA_plot3D.html")
  export_plotly_to_html(pca_plotly_3d_obj, pca_plotly_3d_path)
  cli::cli_alert_success(
    "Interactive 3D PCA saved: {.path {pca_plotly_3d_path}}"
  )

  # ---- 11. Save RDS ----
  cli::cli_h1("Saving session")
  rds <- list(results, comparisons, out_dirs, pca_plot,
              pca_plotly_2d_obj, pca_plotly_3d_obj, annotation,
              padj, fold_change)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  rds_name <- paste0("analysis_results_", ts, ".rds")
  rds_path <- file.path(outdir, rds_name)
  saveRDS(rds, rds_path)
  cli::cli_alert_success("RDS saved: {.path {rds_path}}")

  # ---- 12. Structured run summary (.report.json) ----
  .write_report_json(
    path             = file.path(out_dirs$logs, paste0(ts, "_report.json")),
    runid            = runid,          brs_ticket       = brs_ticket,
    annotation       = annotation,     id_type          = id_type,
    counts_path      = counts_path,    samplesheet_path = samplesheet_path,
    outdir           = outdir,
    padj             = padj,           lfc              = lfc,
    started          = run_started,    finished         = Sys.time(),
    n_genes_input    = nrow(counts),   n_genes_filtered = nrow(dds),
    n_samples        = ncol(countsdf), n_samplesheet_rows = nrow(samplesheet),
    comparisons      = comparisons,    results          = results,
    rds_path         = rds_path,
    session_log      = if (!is.null(session_log)) session_log$path else NA_character_
  )

  cli::cli_h1("Pipeline complete")

  invisible(list(
    results       = results,
    comparisons   = comparisons,
    out_dirs      = out_dirs,
    pca_plot      = pca_plot,
    pca_plotly    = pca_plotly_2d_obj,
    pca_plotly_3d = pca_plotly_3d_obj,
    annotation    = annotation,
    runid         = runid,
    brs_ticket    = brs_ticket,
    sample_info   = sample_info,
    tmm           = tmm,
    counts        = countsdf,
    dds           = dds,
    rds_path      = rds_path
  ))
}
