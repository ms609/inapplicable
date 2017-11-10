library(testthat)

context("Test tree randomness")
test_that("four-tip trees are randomly distributed", {
  nTrees <- 36000
  stringency <- 0.005
  nTip <- 4
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees, 1/(nTip - 1))
  rTrees <- vapply(logical(nTrees), function (XX) unlist(RandomMorphyTree(nTip)), integer((nTip * 4) - 3))
  expect_true(all(rTrees[1 + (seq_len(nTip - 1)), ] %in% nTip + seq_len(nTip - 2)))
  expect_true(expectedBounds[1] < sum(rTrees[2, ] == 5) && expectedBounds[2] > sum(rTrees[2, ] == 5))
  expect_true(expectedBounds[1] < sum(rTrees[3, ] == 5) && expectedBounds[2] > sum(rTrees[3, ] == 5))
  expect_true(expectedBounds[1] < sum(rTrees[4, ] == 5) && expectedBounds[2] > sum(rTrees[4, ] == 5))

  expect_true(all(table(rTrees[c(9, 12), ])[seq_len(nTip - 1)] > expectedBounds[1]))
  expect_true(all(table(rTrees[c(9, 12), ])[seq_len(nTip - 1)] < expectedBounds[2]))

  expect_true(all(table(rTrees[c(10, 13), ])[seq_len(nTip - 1)] < nTrees - expectedBounds[1]))
  expect_true(all(table(rTrees[c(10, 13), ])[seq_len(nTip - 1)] > nTrees - expectedBounds[2]))
})
