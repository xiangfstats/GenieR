context("Phylo input")

test_that("To ensure the valid input", {
  expect_error(branching.sampling.times(matrix(2,2,2)))
})
