library(testthat)

context("Test tree randomness")
test_that("four-tip trees are randomly distributed", {
  nTrees <- 36000
  stringency <- 0.005 # low numbers mean you'll rarely fail by chance
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

test_that("four-tip trees are randomly scored", {
  nTrees <- 36000
  stringency <- 0.005
  nTip <- 4
  
  morphyObj <- mpl_new_Morphy()
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  expect_equal(0, mpl_init_Morphy(nTip, 1, morphyObj) -> error) 
  expect_equal(0, mpl_attach_rawdata('0011;', morphyObj) -> error) 
  expect_equal(0, mpl_set_num_internal_nodes(nTip - 1L, morphyObj) -> error) 
  expect_equal(0, mpl_set_parsim_t(1, 'FITCH', morphyObj) -> error)
  expect_equal(0, mpl_set_charac_weight(1, 1, morphyObj) -> error) 
  expect_equal(0, mpl_apply_tipdata(morphyObj) -> error) 
  class(morphyObj) <- 'morphyPtr'

  scores <- vapply(logical(nTrees), function (XX) {
    testTree <- RandomMorphyTree(nTip)
    MorphyLength(testTree[[1]], testTree[[2]], testTree[[3]], morphyObj)
  }, integer(1))
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees, 1/(nTip - 1))
  rScores <- vapply(logical(nTrees), function (XX) RandomTreeScore(nTip, morphyObj), integer(1))
  expect_true(expectedBounds[1] < sum(scores==1) && expectedBounds[2] > sum(scores==1))
})

test_that("six-tip trees are randomly distributed", {
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
