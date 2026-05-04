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
                         id_type = "ensembl") {
  tryCatch(
    {
      cli::cli_h2("Comparison: {.strong {comparison$name}}")
      cli::cli_inform(c(
        "*" = "Experimental: {.val {comparison$exp}}",
        "*" = "Control:      {.val {comparison$ctrl}}"
      ))

      deseq_results <- perform_deseq2_analysis(dds, comparison$exp, comparison$ctrl)
      if (is.null(deseq_results)) {
        cli::cli_alert_warning("DESeq2 analysis returned NULL; skipping comparison")
        return(NULL)
      }

      cli::cli_alert_info("Annotating results...")
      annotated_results <- annotate_results(
        deseq_results,
        id_type       = id_type,
        annotation_db = annotation_db
      )

      output_file <- create_file_path(out_dirs$de_data, "DESeq2_", comparison$name)
      cli::cli_alert_info("Saving DE results: {.path {output_file}}")
      utils::write.csv(annotated_results, output_file)

      cli::cli_alert_info("Generating volcano plot...")
      volcano_plot <- generate_volcano(
        annotated_results,
        comparison$exp,
        comparison$ctrl
      )
      save_plot(
        volcano_plot,
        create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png")
      )

      cli::cli_alert_info("Generating heatmap (all significant DEGs)...")
      grDevices::png(
        create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap_all_sig.png"),
        width = 800, height = 1200, res = 150
      )
      generate_heatmap(annotated_results, normalized_counts, sample_info,
                       exp_name = comparison$exp, ctrl_name = comparison$ctrl)
      grDevices::dev.off()

      cli::cli_alert_info("Generating heatmap (top 100 DEGs by padj)...")
      grDevices::png(
        create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap_top100.png"),
        width = 800, height = 1400, res = 150
      )
      generate_heatmap(annotated_results, normalized_counts, sample_info,
                       exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                       top_n = 100)
      grDevices::dev.off()

      # ---- Enrichment: GO ----
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
                         id_type    = "ensembl",
                         brs_ticket = "") {
  annotation <- match.arg(annotation, c("human", "mouse"))
  id_type    <- match.arg(id_type, c("ensembl", "entrez", "symbol"))

  # ---- 1. Output dirs + session log ----
  out_dirs <- setup_directories(outdir)
  session_log <- start_session_log(out_dirs$logs)
  on.exit(stop_session_log(session_log), add = TRUE)

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

  # ---- 4. Parse contrasts ----
  cli::cli_h1("Parsing contrasts")
  comparisons <- parse_contrasts(samplesheet)

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

  # ---- 8. Per-comparison loop ----
  cli::cli_h1("Per-comparison differential expression")
  results <- vector("list", length(comparisons))

  cli::cli_progress_bar(
    name   = "Comparisons",
    total  = length(comparisons),
    format = paste0(
      "{cli::pb_name} {cli::pb_current}/{cli::pb_total} | ",
      "{cli::pb_extra$current} | ETA {cli::pb_eta} | ",
      "[{cli::pb_bar}] {cli::pb_percent}"
    ),
    extra  = list(current = ""),
    clear  = FALSE
  )

  for (i in seq_along(comparisons)) {
    cli::cli_progress_update(extra = list(current = comparisons[[i]]$name))
    res <- run_analysis(
      comparison        = comparisons[[i]],
      dds               = dds,
      normalized_counts = tmm,
      sample_info       = sample_info,
      out_dirs          = out_dirs,
      annotation        = annotation,
      annotation_db     = annotation_db,
      id_type           = id_type
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

  sig_genes_union <- unique(unlist(lapply(results, function(r) {
    if (is.null(r) || is.null(r$deseq)) return(character(0))
    d <- r$deseq
    ok <- !is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 0.58
    if (!any(ok)) return(character(0))
    if ("ENSEMBL_ID" %in% colnames(d)) as.character(d$ENSEMBL_ID[ok])
    else rownames(d)[ok]
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
              pca_plotly_2d_obj, pca_plotly_3d_obj, annotation)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  rds_name <- paste0("analysis_results_", ts, ".rds")
  rds_path <- file.path(outdir, rds_name)
  saveRDS(rds, rds_path)
  cli::cli_alert_success("RDS saved: {.path {rds_path}}")

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
