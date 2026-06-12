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
#' @return Annotated data frame that always carries `ENSEMBL_ID`,
#'   `ENTREZID`, `SYMBOL`, and `GENENAME` columns regardless of the input
#'   `id_type`, plus the original DESeq2 result columns.
#' @details Resolves ambiguous one-to-many mappings (e.g. a Symbol that
#'   maps to multiple Ensembl IDs) by keeping the first match and emitting
#'   a `cli::cli_alert_warning`. Pre-checks the input rownames for
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
                                                  envir = .GlobalEnv)) {
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
      "annotate_results: {n_amb} input ID{?s} have ambiguous {.val {id_type}} -> ortholog mappings; keeping first match"
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
  }

  merged
}
