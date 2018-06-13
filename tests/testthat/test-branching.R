context("Phylo input/output")

test_that("To ensure the valid input", {
  expect_error(branching.sampling.times(matrix(2,2,2)))
})

test_that("To ensure the correct output", {
  library(ape)
  t1=rcoal(20)
  x=att(t1)
  expect_true(class(x)=="data.frame" && length(x)==2)
})
