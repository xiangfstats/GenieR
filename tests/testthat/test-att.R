context("Phylo and eps input")

test_that("To ensure the valid input", {
  expect_error(att(matrix(2,2,2)))
})

test_that("effect of changing eps", {
  library(genieR)
  data(village)
  expect_false(dim(att(village))[1]==dim(att(village,eps=1))[1])
  })
