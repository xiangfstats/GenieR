context("Phylo input")

test_that("To ensure the valid input", {
  expect_error(heterochronous.gp.stat(matrix(2,2,2)))
})

test_that("To ensure the correct output", {
  library(ape)
  t1=rcoal(20)
  x=att(t1)
  expect_true(class(heterochronous.gp.stat(t1))=="list" && length(heterochronous.gp.stat(t1))==3)
})
