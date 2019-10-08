context("tree_numbering.R")
test_that("replacement reorder functions work correctly", {
  ## Tree
  tree <- ape::read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  expect_equal(Cladewise(tree), ape::reorder.phylo(tree, 'cladewise'))
  expect_equal(Pruningwise(tree), ape::reorder.phylo(tree, 'pruningwise'))
  expect_equal(Postorder(tree), ape::reorder.phylo(tree, 'postorder'))
})

test_that('ReorderTreeStrict works', {
  tree1 <- ape::read.tree(text = '(((1, 2), 3), (4, (5, 6)));')
  tree2 <- RenumberTips(ape::read.tree(text = '(((6, 5), 4), (3, (1, 2)));'), 
                        as.character(1:6))
  tree3 <- RenumberTips(ape::read.tree(text = '((4, (5, 6)), ((1, 2), 3));'),
                        as.character(1:6))
  
  expect_equal(
    RenumberTreeStrict(tree1$edge[, 1], tree1$edge[, 2]),
    RenumberTreeStrict(tree3$edge[, 1], tree3$edge[, 2]))
  expect_equal(
    RenumberTreeStrict(tree2$edge[, 1], tree2$edge[, 2]),
    RenumberTreeStrict(tree3$edge[, 1], tree3$edge[, 2]))
})