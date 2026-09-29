# The report template is code that no unit test exercised, and v1.7.0 shipped
# with a child template that could not render a run with two or more
# comparisons. Two guards:
#
#   1. A structural check, always on: children under inst/qmd/_sections/ are
#      knit once per comparison, so a labelled chunk collides on the second
#      pass under Quarto ("Duplicate chunk label"). Plain knitr does not
#      reproduce the failure, so the rule is checked on the source text.
#   2. An end-to-end render of the bundled multi-comparison example. It runs
#      the whole pipeline (one to two minutes, KEGG over the network), so it
#      is opt-in: BISR_TEST_RENDER=1, or `just test-render`.
#
# The structural check covers one known cause. Only a render proves that a
# report can be produced, which is what `just release-check` is for.

# Split chunk options on the commas that separate them, leaving alone the
# commas inside strings and calls: fig.cap = "a, b" and fig.dim = c(4, 3).
# A backslash inside a string escapes the next character, so an escaped
# quote does not end the string.
.split_options <- function(x) {
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  quote <- ""
  depth <- 0L
  start <- 1L
  out   <- character(0)
  skip  <- FALSE
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (skip) {
      skip <- FALSE
    } else if (nzchar(quote)) {
      if (ch == "\\") skip <- TRUE else if (ch == quote) quote <- ""
    } else if (ch %in% c("\"", "'")) {
      quote <- ch
    } else if (ch %in% c("(", "[")) {
      depth <- depth + 1L
    } else if (ch %in% c(")", "]")) {
      depth <- depth - 1L
    } else if (ch == "," && depth == 0L) {
      out   <- c(out, substr(x, start, i - 1L))
      start <- i + 1L
    }
  }
  trimws(c(out, substr(x, start, nchar(x))))
}

# Labels declared by the chunks of a knitr/Quarto document, given its lines.
#
# Header: follows knitr's rules, checked against knitr 1.50. The pattern is
# knitr::all_patterns$md$chunk.begin, leading commas are dropped, and ANY
# option without a name is the label wherever it sits, so ```{r setup},
# ```{r, setup} and ```{r, echo=FALSE, setup} all label the chunk, as does
# label="setup".
#
# Option comments under the header (#| label: setup, #| id: setup, the
# comma-separated form, a value on the next line, other engines' prefixes):
# read by knitr::partition_chunk(), the parser knitr itself uses.
.chunk_labels <- function(lines) {
  header <- regmatches(lines, regexec(
    "^[\t >]*```+\\s*\\{([a-zA-Z0-9_]+( *[ ,].*)?)\\}\\s*$", lines))
  at     <- which(lengths(header) > 1)
  opts   <- vapply(header[at], `[[`, character(1), 2)
  engine <- sub("^([a-zA-Z0-9_]+).*$", "\\1", opts)
  opts   <- sub("^[a-zA-Z0-9_]+", "", opts)
  opts   <- gsub("^[[:space:],]+|[[:space:],]+$", "", opts)

  from_header <- unlist(lapply(opts, function(o) {
    tokens <- .split_options(o)
    tokens <- tokens[nzchar(tokens)]
    # an "=" inside a string does not name an option
    named  <- grepl("=", gsub("\"[^\"]*\"|'[^']*'", "", tokens), fixed = TRUE)
    is_label <- grepl("^label\\s*=", tokens)
    c(tokens[!named],
      gsub("^[\"']|[\"']$", "", trimws(sub("^label\\s*=", "", tokens[is_label]))))
  }))

  from_comments <- unlist(lapply(seq_along(at), function(k) {
    body   <- lines[seq_len(length(lines) - at[k]) + at[k]]
    indent <- sub("^([\t >]*).*$", "\\1", lines[at[k]])
    body   <- substring(body, nchar(indent) + 1L)
    knitr::partition_chunk(engine[k], body)$options$label
  }))

  c(as.character(from_header), as.character(from_comments))
}

test_that("the chunk-label detector agrees with knitr on what labels a chunk", {
  labelled <- c(
    "```{r setup}"                                 = "setup",
    "```{r shrink-this-comparison, include=FALSE}" = "shrink-this-comparison",
    "```{r, setup}"                                = "setup",
    "```{r , setup, echo=FALSE}"                   = "setup",
    "```{r, echo=FALSE, setup}"                    = "setup",
    "```{r, label=\"named\", echo=FALSE}"          = "named",
    "```{r, label = named}"                        = "named",
    "> ```{r quoted}"                              = "quoted",
    "\t```{r tabbed}"                              = "tabbed",
    "````{r four}"                                 = "four",
    "``` {r spaced}"                               = "spaced",
    "  ```{python plot}"                           = "plot"
  )
  for (h in names(labelled)) {
    expect_identical(.chunk_labels(h), labelled[[h]], info = h)
  }
  in_comments <- list(
    c("```{r}", "#| label: yaml-label", "1 + 1", "```"),
    c("```{r}", "#| id: yaml-label", "1 + 1", "```"),
    c("```{r}", "#| label=\"yaml-label\", echo=FALSE", "1 + 1", "```"),
    c("```{r}", "#| yaml-label, echo=FALSE", "1 + 1", "```"),
    c("```{r}", "#| label:", "#|   yaml-label", "1 + 1", "```"),
    c("```{r}", "#| echo: false", "#| label: yaml-label", "1 + 1", "```"),
    c("  ```{r}", "  #| label: yaml-label", "  1 + 1", "  ```")
  )
  for (chunk in in_comments) {
    expect_identical(.chunk_labels(chunk), "yaml-label",
                     info = paste(chunk, collapse = " / "))
  }

  unlabelled <- c("```{r}", "```{r, include=FALSE}", "```{r include=FALSE}",
                  "```{r, results='asis', echo=FALSE}",
                  "```{r fig.cap=\"a, b\", echo=FALSE}",
                  "```{r, fig.cap=\"Volcano, shrunken\"}",
                  "```{r, fig.dim=c(4, 3)}",
                  "```{r, fig.cap=\"He said \\\"go\\\", then left\"}",
                  "```", "```r", "~~~{r tilde}",
                  "Inline `r 1 + 1` and a fence ```{r} mid-line are not headers.",
                  "#| echo: false")
  for (h in unlabelled) {
    expect_identical(.chunk_labels(h), character(0), info = h)
  }

  # An option comment counts only at the top of a chunk.
  expect_identical(
    .chunk_labels(c("Prose that mentions", "#| label: not-a-chunk-option", "",
                    "```{r}", "x <- 1", "#| label: too-late", "```")),
    character(0))
})

test_that("per-comparison child templates contain no labelled chunks", {
  sections <- system.file("qmd/_sections", package = "bisrDE")
  expect_true(nzchar(sections) && dir.exists(sections))

  children <- list.files(sections, pattern = "\\.qmd$", full.names = TRUE)
  expect_gt(length(children), 0)

  for (child in children) {
    labels <- .chunk_labels(readLines(child, warn = FALSE))
    expect_identical(
      labels, character(0),
      info = paste0(
        basename(child), " labels a chunk (", paste(labels, collapse = ", "),
        "). report.qmd knits this file once per comparison, so under Quarto a ",
        "label collides on the second comparison and no report is produced."))
  }
})

test_that("the report renders for a run with several comparisons", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("BISR_TEST_RENDER")),
              "opt-in: set BISR_TEST_RENDER=1 (or run `just test-render`)")
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()), "the Quarto CLI is not on PATH")

  assets      <- test_path("..", "..", "..", "assets")
  counts      <- file.path(assets, "example_counts.tsv")
  samplesheet <- file.path(assets, "example_samplesheet.csv")
  skip_if_not(file.exists(counts) && file.exists(samplesheet),
              "the bundled example data is only reachable from a source checkout")

  outdir <- tempfile("bisrde_render_test_")

  # Ten genes make DESeq2 and ggrepel warn (factor conversion, labels dropped
  # outside the plot range). The unit tests own the pipeline's behaviour; this
  # test owns the render, so those warnings are noise here.
  res <- suppressWarnings(suppressMessages(run_pipeline(
    counts_path      = normalizePath(counts),
    samplesheet_path = normalizePath(samplesheet),
    outdir           = outdir,
    runid            = "render_test",
    annotation       = "mouse"
  )))

  # The point of the test is the SECOND pass through the child template. If
  # the example ever shrinks to one comparison the render proves nothing.
  comparisons <- vapply(res$comparisons, function(cmp) cmp$name, character(1))
  expect_gte(length(comparisons), 2)

  # run_pipeline() survives a comparison whose analysis failed, and the
  # template still writes that comparison's heading, so a report full of
  # empty sections would otherwise pass.
  failed <- comparisons[vapply(res$results, function(r) {
    is.null(r) || is.null(r$deseq)
  }, logical(1))]
  expect_identical(failed, character(0),
                   info = "these comparisons produced no DE result")

  html <- suppressMessages(generate_report(rds_path = res$rds_path,
                                           output_dir = outdir))
  expect_true(file.exists(html))

  page     <- readLines(html, warn = FALSE)
  headings <- unlist(regmatches(page, gregexpr("<h3[^>]*>.*?</h3>", page)))
  headings <- trimws(gsub("<[^>]+>", "", headings))
  for (name in comparisons) {
    expect_true(name %in% headings,
                info = paste("no section heading for comparison", name))
  }
})
