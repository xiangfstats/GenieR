context("Phylo input")

test_that("To ensure the valid input", {
  expect_error(heterochronous.gp.stat(matrix(2,2,2)))
})
