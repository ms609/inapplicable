library(ape)

context("tree search")
test_that("tree search finds shortest tree", {
  ## Tree
  tree <- read.tree(text = "(((((1,2),3),4),5),6);")
  dataset <- StringToPhyDat('110000 111000 111100', 1:6, byTaxon=FALSE)
  expect_equal(Cladewise(tree), reorder(tree, 'cladewise'))
  expect_equal(Pruningwise(tree), reorder(tree, 'pruningwise'))
  expect_equal(Postorder(tree), reorder(tree, 'postorder'))
})
