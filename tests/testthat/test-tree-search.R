library(ape)

context("tree search")
test_that("tree search finds shortest tree", {
  true_tree <- read.tree(text = "(((((1,2),3),4),5),6);")
  dataset <- StringToPhyDat('110000 111000 111100', 1:6, byTaxon=FALSE)
  start_tree <- TreeSearch::RenumberTips(read.tree(text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)
  expect_equal(InapplicableFitch(start_tree, dataset), 6)
  expect_equal(InapplicableFitch(TreeSearch(start_tree, dataset, method='NNI', verbosity=0), dataset), InapplicableFitch(true_tree, dataset), 3)
  expect_equal(InapplicableFitch(TreeSearch(start_tree, dataset, method='SPR', verbosity=0), dataset), InapplicableFitch(true_tree, dataset), 3)
  expect_equal(InapplicableFitch(TreeSearch(start_tree, dataset, method='TBR', verbosity=0), dataset), InapplicableFitch(true_tree, dataset), 3)
  expect_equal(InapplicableFitch(Ratchet(start_tree, dataset, k=3, maxhits=5, verbosity=0), dataset), InapplicableFitch(true_tree, dataset), 3)
})
