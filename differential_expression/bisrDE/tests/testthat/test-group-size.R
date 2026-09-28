# v1.7: the pre-filter's smallest_group_size is derived from the design.

test_that(".smallest_group_size follows the smallest group with a floor of 2", {
  si <- data.frame(sample = paste0("S", 1:9),
                   condition = c(rep("A", 4), rep("B", 3), rep("C", 2)))
  expect_equal(bisrDE:::.smallest_group_size(si), 2L)
  si2 <- data.frame(sample = paste0("S", 1:8), condition = rep(c("A", "B"), each = 4))
  expect_equal(bisrDE:::.smallest_group_size(si2), 4L)
  si3 <- data.frame(sample = c("S1", "S2"), condition = c("A", "B"))
  expect_equal(bisrDE:::.smallest_group_size(si3), 2L)   # floor
})
