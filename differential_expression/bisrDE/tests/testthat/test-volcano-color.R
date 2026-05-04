# Phase 1.1 fix verification: the volcano plot's color tags, threshold
# line, and y-axis must all use the same significance field (`sig`),
# and the 3 color levels must map to the documented palette (Over =
# firebrick, Under = steelblue, Not significant = grey70).

build_volcano_data <- function() {
  data.frame(
    ENSEMBL_ID     = paste0("ENSG", sprintf("%07d", 1:6)),
    SYMBOL         = paste0("Gene", 1:6),
    log2FoldChange = c(2.0, -2.0, 0.0, 1.0, -1.0, 0.1),
    pvalue         = c(1e-4, 1e-4, 0.5, 1e-4, 1e-4, 1e-4),
    padj           = c(1e-3, 1e-3, 0.5, 1e-3, 1e-3, 1e-3),
    stringsAsFactors = FALSE
  )
}

test_that("generate_volcano returns a ggplot with the 3-level color_tag factor", {
  p <- generate_volcano(build_volcano_data(),
                        exp_name = "A", ctrl_name = "B")

  expect_s3_class(p, "ggplot")
  expect_equal(
    levels(p$data$color_tag),
    c("Over expressed", "Under expressed", "Not significant")
  )
})

test_that("generate_volcano color_tag classifies up/down/non-sig correctly", {
  p <- generate_volcano(build_volcano_data(),
                        exp_name = "A", ctrl_name = "B",
                        p = 0.05, lfc = 0.58)

  # Per the fixture:
  #  row 1: lfc =  2.0,  padj = 1e-3 -> "Over expressed"
  #  row 2: lfc = -2.0,  padj = 1e-3 -> "Under expressed"
  #  row 3: lfc =  0.0,  padj = 0.5  -> "Not significant"
  #  row 4: lfc =  1.0,  padj = 1e-3 -> "Over expressed"
  #  row 5: lfc = -1.0,  padj = 1e-3 -> "Under expressed"
  #  row 6: lfc =  0.1,  padj = 1e-3 -> "Not significant" (|lfc| < 0.58)
  expect_equal(
    as.character(p$data$color_tag),
    c("Over expressed", "Under expressed", "Not significant",
      "Over expressed", "Under expressed", "Not significant")
  )
})

test_that("generate_volcano uses the documented palette", {
  p <- generate_volcano(build_volcano_data(),
                        exp_name = "A", ctrl_name = "B")

  # Find the colour scale (scale_color_manual)
  color_scale <- NULL
  for (s in p$scales$scales) {
    if ("colour" %in% s$aesthetics) {
      color_scale <- s
      break
    }
  }
  expect_false(is.null(color_scale))

  # The palette stores the named-vector `values =` argument we passed
  expect_equal(color_scale$palette(3),
               c("Over expressed"  = "firebrick",
                 "Under expressed" = "steelblue",
                 "Not significant" = "grey70"))
})

test_that("generate_volcano horizontal threshold line is at -log10(p)", {
  p <- generate_volcano(build_volcano_data(),
                        exp_name = "A", ctrl_name = "B",
                        p = 0.05, lfc = 0.58)

  # Find the geom_hline layer (its data has yintercept).
  hline_layer <- NULL
  for (l in p$layers) {
    if (inherits(l$geom, "GeomHline")) {
      hline_layer <- l
      break
    }
  }
  expect_false(is.null(hline_layer))
  yint <- hline_layer$data$yintercept[1]
  expect_equal(yint, -log10(0.05))
})

test_that("generate_volcano vertical threshold lines are at +/- lfc", {
  p <- generate_volcano(build_volcano_data(),
                        exp_name = "A", ctrl_name = "B",
                        p = 0.05, lfc = 0.58)

  vline_layer <- NULL
  for (l in p$layers) {
    if (inherits(l$geom, "GeomVline")) {
      vline_layer <- l
      break
    }
  }
  expect_false(is.null(vline_layer))
  expect_setequal(vline_layer$data$xintercept, c(-0.58, 0.58))
})

test_that("generate_volcano y-axis labels use the chosen sig column", {
  p <- generate_volcano(build_volcano_data(),
                        exp_name = "A", ctrl_name = "B",
                        sig = "padj")

  expect_equal(p$labels$y, "-log10( padj )")

  p2 <- generate_volcano(build_volcano_data(),
                         exp_name = "A", ctrl_name = "B",
                         sig = "pvalue")

  expect_equal(p2$labels$y, "-log10( pvalue )")
})

test_that("generate_volcano title shows 'exp vs. ctrl'", {
  p <- generate_volcano(build_volcano_data(),
                        exp_name = "Treated", ctrl_name = "Control")
  expect_match(p$labels$title, "Treated vs\\. Control")
})
