context("Phylo and eps input/output")

test_that("To ensure the valid input", {
  expect_error(att(matrix(2,2,2)))
})

test_that("effect of changing eps", {
  library(genieR)
  data(village)
  expect_false(dim(att(village))[1]==dim(att(village,eps=1))[1])
  })

test_that("To ensure the correct output", {
  library(ape)
  t1=rcoal(20)
  x=att(t1)
  expect_true(dim(x)[2]==2)
})
