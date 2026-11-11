# Samplesheet header handling. A `Sample_ID,Group_ID` header used to parse
# "fine" but leave samplesheet$GroupID NULL, so every contrast silently found
# no groups and the run died much later with the misleading
# "No contrasts to analyse after include/exclude filtering and group checks".

write_ss <- function(lines) {
  p <- tempfile(fileext = ".csv")
  writeLines(lines, p)
  p
}

test_that("read_samplesheet normalizes Sample_ID / Group_ID header variants", {
  p <- write_ss(c(
    "Sample_ID,Group_ID,OE_RSF1_vs_GFP",
    "s1,GFP,0",
    "s2,OE_RSF1,1"
  ))

  ss <- read_samplesheet(p)

  expect_true(all(c("SampleID", "GroupID") %in% colnames(ss)))
  expect_false(any(c("Sample_ID", "Group_ID") %in% colnames(ss)))
  # the contrast column is untouched
  expect_true("OE_RSF1_vs_GFP" %in% colnames(ss))
})

test_that("a normalized samplesheet actually yields comparisons", {
  # End-to-end regression for the reported failure.
  p <- write_ss(c(
    "Sample_ID,Group_ID,OE_RSF1_vs_GFP",
    "s1,GFP,0",
    "s2,GFP,0",
    "s3,OE_RSF1,1",
    "s4,OE_RSF1,1"
  ))

  cmp <- parse_contrasts(read_samplesheet(p))

  expect_length(cmp, 1)
  expect_equal(cmp[[1]]$name, "OE_RSF1_vs_GFP")
  expect_equal(cmp[[1]]$exp,  "OE_RSF1")
  expect_equal(cmp[[1]]$ctrl, "GFP")
})

test_that("read_samplesheet accepts other common variants", {
  p <- write_ss(c("sample id,condition,a_vs_b", "s1,A,1", "s2,B,0"))
  ss <- read_samplesheet(p)
  expect_true(all(c("SampleID", "GroupID") %in% colnames(ss)))
})

test_that("canonical headers are left exactly as-is", {
  p <- write_ss(c("SampleID,GroupID,a_vs_b", "s1,A,1", "s2,B,0"))
  ss <- read_samplesheet(p)
  expect_equal(colnames(ss), c("SampleID", "GroupID", "a_vs_b"))
})

test_that("read_samplesheet errors clearly when GroupID is absent entirely", {
  p <- write_ss(c("SampleID,Batch,a_vs_b", "s1,1,1", "s2,2,0"))
  expect_error(read_samplesheet(p), "GroupID")
})

test_that("normalization never renames a contrast column or duplicates a name", {
  # A contrast column whose name contains an alias token must survive intact...
  p <- write_ss(c("SampleID,GroupID,group_vs_control", "s1,A,1", "s2,B,0"))
  ss <- read_samplesheet(p)
  expect_true("group_vs_control" %in% colnames(ss))
  expect_equal(sum(colnames(ss) == "GroupID"), 1L)

  # ...and a sheet carrying BOTH spellings must not gain a duplicate column.
  p2 <- write_ss(c("SampleID,Sample_ID,GroupID,a_vs_b", "s1,x,A,1", "s2,y,B,0"))
  ss2 <- read_samplesheet(p2)
  expect_equal(sum(colnames(ss2) == "SampleID"), 1L)
})
