# Annotation: translate AnnotationDbi keytypes and merge gene-ID columns
# onto DESeq2 result tables. Source-of-truth during Phase 5.2-5.7 is still
# the parent's de.R; this package copy is a snapshot per the Phase 5.3 spec.

#' Translate a `--id-type` value to an AnnotationDbi keytype
#'
#' @description Internal switch from the user-facing `id_type` strings
#'   ("ensembl", "entrez", "symbol") to the AnnotationDbi keytype names
#'   ("ENSEMBL", "ENTREZID", "SYMBOL").
#'
#' @param id_type One of "ensembl", "entrez", "symbol" (case-insensitive).
#' @return The corresponding AnnotationDbi keytype string.
#' @keywords internal
.id_type_to_keytype <- function(id_type) {
  switch(tolower(id_type),
    "ensembl" = "ENSEMBL",
    "entrez"  = "ENTREZID",
    "symbol"  = "SYMBOL",
    stop("Invalid id_type: '", id_type, "' (must be 'ensembl', 'entrez', or 'symbol')")
  )
}


#' Annotate DESeq2 results with gene identifier columns
#'
#' @description Join Entrez, Symbol, GENENAME, and Ensembl ID columns onto
#'   a DESeq2 results data frame using the org.*.eg.db package matching
#'   `--annotation`. The output always carries an `ENSEMBL_ID` column
#'   regardless of the input `id_type`, so downstream plotting / report
#'   code that hardcodes `ENSEMBL_ID` keeps working.
#'
#' @param results DESeq2 results data frame whose rownames are gene
#'   identifiers of type `id_type` (typically Ensembl gene IDs).
#' @param id_type Identifier type used in `rownames(results)`. Defaults to
#'   the `id_type` global set during option parsing in the parent
#'   pipeline; falls back to `"ensembl"` if no global is set. One of
#'   "ensembl", "entrez", "symbol".
#' @param annotation_db An `OrgDb` object (e.g. `org.Hs.eg.db::org.Hs.eg.db`
#'   or `org.Mm.eg.db::org.Mm.eg.db`). Defaults to the `annotation_db`
#'   global set during option parsing in the parent pipeline. Errors if
#'   neither the argument nor the global is supplied.
#' @param gene_meta Optional data frame with columns `id`, `id_versioned`,
#'   `gene_name` (the `"gene_meta"` attribute of [read_counts()]). When
#'   supplied: genes the OrgDb cannot name get their `SYMBOL` from the
#'   input's `gene_name` (recorded in `SYMBOL_SOURCE` as `"input"` vs
#'   `"orgdb"`), and an `ENSEMBL_ID_VERSIONED` column carries the original
#'   versioned accession when any input ID had a version suffix. Default
#'   `NULL` (no fallback, no versioned column).
#' @return Annotated data frame that always carries `ENSEMBL_ID`,
#'   `ENTREZID`, `SYMBOL`, and `GENENAME` columns regardless of the input
#'   `id_type`, plus the original DESeq2 result columns (and, with
#'   `gene_meta`, `SYMBOL_SOURCE` / `ENSEMBL_ID_VERSIONED`).
#' @details One-to-many mappings (one input ID with several OrgDb rows, e.g.
#'   an Ensembl gene with two Entrez IDs, or a Symbol matching several
#'   Ensembl genes) are collapsed by keeping the FIRST row `AnnotationDbi`
#'   returns; a `cli::cli_alert_warning` reports how many. This is a
#'   deterministic but arbitrary rule and is stated in the report's Methods.
#'   Pre-checks the input rownames for
#'   duplicates and emits a `cli::cli_alert_danger` if found (treated as a
#'   data-quality problem rather than a hard stop, so the pipeline can
#'   still finish - but the user is alerted). When the input `id_type` is
#'   not Ensembl, the original input ID is preserved under
#'   `{KEYTYPE}_ID` and the AnnotationDbi-mapped Ensembl value is promoted
#'   to `ENSEMBL_ID` for downstream-code compatibility.
#'
#' @importFrom dplyr rename
#' @export
annotate_results <- function(results,
                             id_type = get0("id_type", envir = .GlobalEnv,
                                            ifnotfound = "ensembl"),
                             annotation_db = get0("annotation_db",
                                                  envir = .GlobalEnv),
                             gene_meta = NULL) {
  if (is.null(annotation_db)) {
    stop("annotate_results: `annotation_db` is NULL. Either pass it ",
         "explicitly or set the `annotation_db` global before calling.")
  }

  keytype <- .id_type_to_keytype(id_type)
  gene_ids <- rownames(results)

  dup_input <- duplicated(gene_ids)
  if (any(dup_input)) {
    n_dup <- sum(dup_input)
    cli::cli_alert_danger(
      "annotate_results: input has {n_dup} duplicate {.val {id_type}} ID{?s}; downstream merges may double-count rows"
    )
  }

  target_cols <- c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME")
  query_cols  <- setdiff(target_cols, keytype)

  annotations <- AnnotationDbi::select(
    annotation_db,
    keys    = unique(gene_ids),
    keytype = keytype,
    columns = query_cols
  )

  dup_mask <- duplicated(annotations[[keytype]])
  if (any(dup_mask)) {
    n_amb <- length(unique(annotations[[keytype]][dup_mask]))
    cli::cli_alert_warning(
      "annotate_results: {n_amb} input ID{?s} map to more than one annotation row (one-to-many {.val {id_type}} mapping); keeping the first match"
    )
    annotations <- annotations[!dup_mask, ]
  }

  merged <- merge(
    as.data.frame(results),
    annotations,
    by.x  = "row.names",
    by.y  = keytype,
    all.x = TRUE
  )

  if (keytype == "ENSEMBL") {
    merged <- merged |> dplyr::rename(ENSEMBL_ID = "Row.names")
  } else {
    input_col <- paste0(keytype, "_ID")
    merged <- merged |> dplyr::rename(!!input_col := "Row.names")
    merged$ENSEMBL_ID <- merged$ENSEMBL
    merged$ENSEMBL <- NULL
    # The input keytype was the merge key, so its canonical column (e.g. SYMBOL
    # for id_type="symbol", ENTREZID for "entrez") is NOT among the queried
    # columns — only the "<keytype>_ID" input copy is. The input IDs *are* that
    # column's values, so mirror them in: the documented contract (and the
    # volcano / heatmap / enrichment code) require ENSEMBL_ID, SYMBOL, ENTREZID,
    # and GENENAME to always be present regardless of id_type.
    if (is.null(merged[[keytype]])) merged[[keytype]] <- merged[[input_col]]
  }

  # ---- input gene_name fallback + versioned ID (from read_counts()) ----
  if (!is.null(gene_meta) && all(c("id", "gene_name") %in% colnames(gene_meta))) {
    key_col <- if (keytype == "ENSEMBL") "ENSEMBL_ID" else paste0(keytype, "_ID")
    idx <- match(as.character(merged[[key_col]]), as.character(gene_meta$id))
    input_name <- as.character(gene_meta$gene_name)[idx]
    ver_name   <- if ("id_versioned" %in% colnames(gene_meta)) as.character(gene_meta$id_versioned)[idx] else NA_character_
    # For symbol-keyed input SYMBOL was mirrored from the input IDs above, so
    # "did the OrgDb know this gene" has to be read off a column the OrgDb
    # actually filled in.
    has_orgdb <- if (keytype == "SYMBOL") {
      (!is.na(merged$ENTREZID) & nzchar(as.character(merged$ENTREZID))) |
        (!is.na(merged$GENENAME) & nzchar(as.character(merged$GENENAME)))
    } else {
      !is.na(merged$SYMBOL) & nzchar(as.character(merged$SYMBOL))
    }
    # nf-core writes gene_name = gene_id when the GTF has no name; an
    # accession is not a name, so it must not become the SYMBOL.
    has_input  <- !is.na(input_name) & nzchar(input_name) &
      input_name != as.character(merged[[key_col]]) &
      (is.na(ver_name) | input_name != ver_name) &
      !grepl("^ENS[A-Z]*G[0-9]+", input_name)
    # Symbol-keyed input: the symbol IS the input, so anything the OrgDb did
    # not recognise is by definition input-sourced.
    merged$SYMBOL_SOURCE <- if (keytype == "SYMBOL") {
      ifelse(has_orgdb, "orgdb", "input")
    } else {
      ifelse(has_orgdb, "orgdb", ifelse(has_input, "input", NA_character_))
    }
    fill <- keytype != "SYMBOL" & !has_orgdb & has_input
    if (any(fill)) {
      merged$SYMBOL[fill] <- input_name[fill]
      cli::cli_alert_info(
        "annotate_results: {sum(fill)} gene{?s} without an OrgDb symbol named from the input gene_name column (see SYMBOL_SOURCE)"
      )
    }
    if (keytype == "ENSEMBL" &&
        any(!is.na(ver_name) & ver_name != as.character(merged[[key_col]]))) {
      merged$ENSEMBL_ID_VERSIONED <- ver_name
    }
  }

  merged
}
