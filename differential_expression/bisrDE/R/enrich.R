# Gene set enrichment analysis (GSEA): GO, KEGG, Reactome, MSigDB Hallmark,
# plus a shared dotplot visualizer. Source-of-truth during Phase 5.2-5.7
# is still the parent's de.R; this package copy is a snapshot per the
# Phase 5.3 spec.

#' Map the --annotation flag value to a KEGG organism code
#'
#' @description Internal switch from the user-facing `annotation` strings
#'   ("human"/"mouse") to the KEGG organism codes ("hsa"/"mmu") expected
#'   by `clusterProfiler::gseKEGG`.
#'
#' @param annotation One of "human" or "mouse" (the values accepted by
#'   `--annotation`).
#' @return KEGG organism code: "hsa" for human, "mmu" for mouse.
#' @keywords internal
.kegg_organism <- function(annotation) {
  switch(tolower(annotation),
    "human" = "hsa",
    "mouse" = "mmu",
    stop("Unsupported --annotation '", annotation, "' for KEGG (need 'human' or 'mouse')")
  )
}


#' Map the --annotation flag value to an msigdbr species string
#'
#' @description Internal switch from "human"/"mouse" to the binomial species
#'   names that `msigdbr::msigdbr` expects.
#'
#' @param annotation One of "human" or "mouse".
#' @return msigdbr species string ("Homo sapiens" or "Mus musculus").
#' @keywords internal
.msigdbr_species <- function(annotation) {
  switch(tolower(annotation),
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus",
    stop("Unsupported --annotation '", annotation, "' for MSigDB (need 'human' or 'mouse')")
  )
}


#' GO Gene Set Enrichment Analysis (all ontologies)
#'
#' @description GSEA over GO BP/MF/CC, ranked by `log2FoldChange`. Wraps
#'   `clusterProfiler::gseGO` with `keyType = "ENSEMBL"`. Returns a
#'   `setReadable`'d result so downstream tables/dotplots show gene
#'   symbols rather than raw Ensembl IDs.
#'
#' @param annotated_result Annotated DESeq2 results frame (output of
#'   [annotate_results()]). Must contain `ENSEMBL_ID`, `log2FoldChange`,
#'   and `baseMean` columns.
#' @param p Adjusted p-value cutoff. Default `1` keeps every result; the
#'   report-side filter applies the user-facing threshold.
#' @param annotation_db `OrgDb` object. Defaults to the `annotation_db`
#'   global set during option parsing in the parent pipeline; errors if
#'   neither the arg nor the global is supplied.
#' @return setReadable'd `gseGO` object, or `NULL` if no enriched GO
#'   terms (or any error during enrichment).
#' @details Filters input to genes with `baseMean > 0`, non-empty
#'   `ENSEMBL_ID`, and non-NA `log2FoldChange` before ranking. Duplicate
#'   Ensembl IDs in the rank vector are resolved by first-occurrence to
#'   keep ordering stable.
#'
#' @importFrom clusterProfiler gseGO setReadable
#' @importFrom dplyr filter
#' @export
process_gsea <- function(annotated_result, p = 1,
                         annotation_db = get0("annotation_db",
                                              envir = .GlobalEnv)) {
  if (is.null(annotation_db)) {
    stop("process_gsea: `annotation_db` is NULL. Either pass it explicitly ",
         "or set the `annotation_db` global before calling.")
  }

  tryCatch(
    {
      set.seed(1000)

      df <- annotated_result |>
        dplyr::filter(.data$baseMean > 0,
                      !is.na(.data$ENSEMBL_ID), nzchar(.data$ENSEMBL_ID),
                      !is.na(.data$log2FoldChange))

      if (nrow(df) == 0) {
        cli::cli_alert_warning("GSEA-GO: no genes with valid ENSEMBL_ID + log2FC; skipping")
        return(NULL)
      }

      gene_list <- df$log2FoldChange
      names(gene_list) <- df$ENSEMBL_ID
      gene_list <- gene_list[!duplicated(names(gene_list))]
      gene_list <- sort(gene_list, decreasing = TRUE)

      gse <- clusterProfiler::gseGO(
        geneList      = gene_list,
        ont           = "ALL",
        minGSSize     = 10,
        maxGSSize     = 1000,
        keyType       = "ENSEMBL",
        pvalueCutoff  = p,
        pAdjustMethod = "fdr",
        OrgDb         = annotation_db
      )

      if (is.null(gse) || nrow(gse@result) == 0) {
        cli::cli_alert_warning("GSEA-GO: no enriched terms found")
        return(NULL)
      }

      clusterProfiler::setReadable(gse, OrgDb = annotation_db, keyType = "ENSEMBL")
    },
    error = function(e) {
      cli::cli_alert_danger("GSEA-GO error: {e$message}")
      NULL
    }
  )
}


#' KEGG pathway Gene Set Enrichment Analysis
#'
#' @description GSEA over KEGG pathways, ranked by `log2FoldChange`. Wraps
#'   `clusterProfiler::gseKEGG` with `keyType = "ncbi-geneid"` (Entrez).
#'   Returns a `setReadable`'d result with gene Symbols.
#'
#' @param annotated_result Annotated DESeq2 results frame (output of
#'   [annotate_results()]). Must contain `ENTREZID`, `log2FoldChange`,
#'   and `baseMean` columns.
#' @param p Adjusted p-value cutoff. Default `1`.
#' @param annotation One of "human" or "mouse"; selects the KEGG organism
#'   code via [.kegg_organism()]. Defaults to the `annotation` global set
#'   during option parsing.
#' @param annotation_db `OrgDb` object for `setReadable` Symbol mapping.
#'   Defaults to the `annotation_db` global set during option parsing.
#' @return setReadable'd `gseKEGG` object, or `NULL` if no enriched
#'   pathways or any error.
#' @details Rank-vector construction matches [process_gsea()]: filter to
#'   non-NA Entrez + log2FC, dedupe by first-occurrence, sort descending.
#'
#' @importFrom clusterProfiler gseKEGG setReadable
#' @importFrom dplyr filter
#' @export
process_kegg_gsea <- function(annotated_result, p = 1,
                              annotation = get0("annotation",
                                                envir = .GlobalEnv),
                              annotation_db = get0("annotation_db",
                                                   envir = .GlobalEnv)) {
  if (is.null(annotation)) {
    stop("process_kegg_gsea: `annotation` is NULL. Either pass it ",
         "explicitly or set the `annotation` global before calling.")
  }
  if (is.null(annotation_db)) {
    stop("process_kegg_gsea: `annotation_db` is NULL. Either pass it ",
         "explicitly or set the `annotation_db` global before calling.")
  }

  tryCatch(
    {
      set.seed(1000)

      df <- annotated_result |>
        dplyr::filter(.data$baseMean > 0,
                      !is.na(.data$ENTREZID), nzchar(.data$ENTREZID),
                      !is.na(.data$log2FoldChange))

      if (nrow(df) == 0) {
        cli::cli_alert_warning("KEGG GSEA: no genes with valid ENTREZID + log2FC; skipping")
        return(NULL)
      }

      gene_list <- df$log2FoldChange
      names(gene_list) <- df$ENTREZID
      gene_list <- gene_list[!duplicated(names(gene_list))]
      gene_list <- sort(gene_list, decreasing = TRUE)

      org <- .kegg_organism(annotation)

      gse <- clusterProfiler::gseKEGG(
        geneList      = gene_list,
        organism      = org,
        keyType       = "ncbi-geneid",
        minGSSize     = 10,
        maxGSSize     = 1000,
        pvalueCutoff  = p,
        pAdjustMethod = "fdr",
        verbose       = FALSE
      )

      if (is.null(gse) || nrow(gse@result) == 0) {
        cli::cli_alert_warning("KEGG GSEA: no enriched pathways found")
        return(NULL)
      }

      clusterProfiler::setReadable(gse, OrgDb = annotation_db, keyType = "ENTREZID")
    },
    error = function(e) {
      cli::cli_alert_danger("KEGG GSEA error: {e$message}")
      NULL
    }
  )
}


#' Reactome pathway Gene Set Enrichment Analysis
#'
#' @description GSEA over Reactome pathways. Wraps `ReactomePA::gsePathway`
#'   directly with the `--annotation` value (Reactome accepts "human"
#'   and "mouse" natively, so no organism-code translation is needed —
#'   unlike KEGG).
#'
#' @param annotated_result Annotated DESeq2 results frame (output of
#'   [annotate_results()]). Must contain `ENTREZID`, `log2FoldChange`,
#'   and `baseMean` columns.
#' @param p Adjusted p-value cutoff. Default `1`.
#' @param annotation One of "human" or "mouse". Defaults to the
#'   `annotation` global set during option parsing.
#' @param annotation_db `OrgDb` for `setReadable` Symbol mapping.
#'   Defaults to the `annotation_db` global.
#' @return setReadable'd `gsePathway` object, or `NULL` if no enriched
#'   pathways or any error.
#'
#' @importFrom ReactomePA gsePathway
#' @importFrom clusterProfiler setReadable
#' @importFrom dplyr filter
#' @export
process_reactome_gsea <- function(annotated_result, p = 1,
                                  annotation = get0("annotation",
                                                    envir = .GlobalEnv),
                                  annotation_db = get0("annotation_db",
                                                       envir = .GlobalEnv)) {
  if (is.null(annotation)) {
    stop("process_reactome_gsea: `annotation` is NULL. Either pass it ",
         "explicitly or set the `annotation` global before calling.")
  }
  if (is.null(annotation_db)) {
    stop("process_reactome_gsea: `annotation_db` is NULL. Either pass it ",
         "explicitly or set the `annotation_db` global before calling.")
  }

  tryCatch(
    {
      set.seed(1000)

      df <- annotated_result |>
        dplyr::filter(.data$baseMean > 0,
                      !is.na(.data$ENTREZID), nzchar(.data$ENTREZID),
                      !is.na(.data$log2FoldChange))

      if (nrow(df) == 0) {
        cli::cli_alert_warning("Reactome GSEA: no genes with valid ENTREZID + log2FC; skipping")
        return(NULL)
      }

      gene_list <- df$log2FoldChange
      names(gene_list) <- df$ENTREZID
      gene_list <- gene_list[!duplicated(names(gene_list))]
      gene_list <- sort(gene_list, decreasing = TRUE)

      gse <- ReactomePA::gsePathway(
        geneList      = gene_list,
        organism      = tolower(annotation),
        minGSSize     = 10,
        maxGSSize     = 1000,
        pvalueCutoff  = p,
        pAdjustMethod = "fdr",
        verbose       = FALSE
      )

      if (is.null(gse) || nrow(gse@result) == 0) {
        cli::cli_alert_warning("Reactome GSEA: no enriched pathways found")
        return(NULL)
      }

      clusterProfiler::setReadable(gse, OrgDb = annotation_db, keyType = "ENTREZID")
    },
    error = function(e) {
      cli::cli_alert_danger("Reactome GSEA error: {e$message}")
      NULL
    }
  )
}


#' MSigDB Hallmark gene set GSEA
#'
#' @description GSEA over the MSigDB Hallmark collection (50 well-curated
#'   gene sets). MSigDB has no dedicated GSEA wrapper the way GO/KEGG/
#'   Reactome do; this builds a `TERM2GENE` table from
#'   `msigdbr::msigdbr(species, collection = "H")` and passes it to
#'   `clusterProfiler::GSEA()`.
#'
#' @param annotated_result Annotated DESeq2 results frame (output of
#'   [annotate_results()]). Must contain `ENTREZID`, `log2FoldChange`,
#'   `baseMean`.
#' @param p Adjusted p-value cutoff. Default `1`.
#' @param annotation One of "human" or "mouse". Defaults to the
#'   `annotation` global.
#' @param annotation_db `OrgDb` for `setReadable` Symbol mapping.
#'   Defaults to the `annotation_db` global.
#' @return setReadable'd GSEA object, or `NULL` if no enriched gene sets
#'   or any error.
#' @details Targets msigdbr ≥10.0 (renv lock has v26.1.0) — uses
#'   `species`/`collection` args and the `ncbi_gene` Entrez column. The
#'   pre-v10 args (`category` + `entrez_gene` column) are deprecated but
#'   still functional; the modern names avoid deprecation warnings on
#'   every run.
#'
#' @importFrom clusterProfiler GSEA setReadable
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr filter
#' @export
process_msigdb_hallmark <- function(annotated_result, p = 1,
                                    annotation = get0("annotation",
                                                      envir = .GlobalEnv),
                                    annotation_db = get0("annotation_db",
                                                         envir = .GlobalEnv)) {
  if (is.null(annotation)) {
    stop("process_msigdb_hallmark: `annotation` is NULL. Either pass it ",
         "explicitly or set the `annotation` global before calling.")
  }
  if (is.null(annotation_db)) {
    stop("process_msigdb_hallmark: `annotation_db` is NULL. Either pass it ",
         "explicitly or set the `annotation_db` global before calling.")
  }

  tryCatch(
    {
      set.seed(1000)

      df <- annotated_result |>
        dplyr::filter(.data$baseMean > 0,
                      !is.na(.data$ENTREZID), nzchar(.data$ENTREZID),
                      !is.na(.data$log2FoldChange))

      if (nrow(df) == 0) {
        cli::cli_alert_warning("MSigDB Hallmark: no genes with valid ENTREZID + log2FC; skipping")
        return(NULL)
      }

      gene_list <- df$log2FoldChange
      names(gene_list) <- df$ENTREZID
      gene_list <- gene_list[!duplicated(names(gene_list))]
      gene_list <- sort(gene_list, decreasing = TRUE)

      species <- .msigdbr_species(annotation)
      msig_h  <- msigdbr::msigdbr(species = species, collection = "H")
      term2gene <- data.frame(
        gs_name = msig_h$gs_name,
        gene    = as.character(msig_h$ncbi_gene),
        stringsAsFactors = FALSE
      )

      gse <- clusterProfiler::GSEA(
        geneList      = gene_list,
        TERM2GENE     = term2gene,
        minGSSize     = 10,
        maxGSSize     = 1000,
        pvalueCutoff  = p,
        pAdjustMethod = "fdr",
        verbose       = FALSE
      )

      if (is.null(gse) || nrow(gse@result) == 0) {
        cli::cli_alert_warning("MSigDB Hallmark: no enriched gene sets found")
        return(NULL)
      }

      clusterProfiler::setReadable(gse, OrgDb = annotation_db, keyType = "ENTREZID")
    },
    error = function(e) {
      cli::cli_alert_danger("MSigDB Hallmark error: {e$message}")
      # msigdbr >= 24 fetches the gene-set archive from Zenodo on first use and
      # caches it. On an HPC compute node with no outbound internet that surfaces
      # as an opaque timeout, so name the actual fix rather than leaving the user
      # to guess. Only for network-ish failures — a genuine analysis error should
      # not be buried under installation advice.
      if (grepl("timeout|timed out|resolve host|connection|curl|download|internet",
                e$message, ignore.case = TRUE)) {
        cli::cli_alert_info(c(
          "This looks like {.pkg msigdbr} failing to download its gene sets.",
          "i" = "Warm the cache once on a machine WITH internet (e.g. a login node):",
          "*" = '{.code Rscript -e \'msigdbr::msigdbr(species = "Homo sapiens", collection = "H")\'}',
          "i" = "The cache lives in {.path {tools::R_user_dir(\"msigdbr\", \"cache\")}} and is reused offline.",
          "i" = "Override the location with the {.envvar R_USER_CACHE_DIR} environment variable."
        ))
      }
      NULL
    }
  )
}


#' Dotplot of GSEA results, split by direction of change
#'
#' @description Visualize a GSEA result object as a top-15 dotplot, split
#'   into activated/suppressed panels by the `.sign` column.
#'
#' @param gse A GSEA result object (e.g. output of [process_gsea()],
#'   [process_kegg_gsea()], [process_reactome_gsea()], or
#'   [process_msigdb_hallmark()]).
#' @param title Plot title.
#' @return A `ggplot` object: dot size = gene count, color =
#'   significance, panels = direction of enrichment.
#'
#' @importFrom enrichplot dotplot
#' @importFrom ggplot2 facet_grid
#' @export
create_dotplot <- function(gse, title) {
  enrichplot::dotplot(
    gse,
    showCategory = 15,
    title        = title,
    split        = ".sign",
    orderBy      = "p.adjust",
    label_format = 31,
    font.size    = 9
  ) +
    ggplot2::facet_grid(. ~ .sign) +
    ggplot2::scale_color_viridis_c(option = "C", direction = -1)
}
