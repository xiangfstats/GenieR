context("loglikelihood of two methods")

test_that("same loglikelihood", {
  library(ape)
  t1=rcoal(20)
  expect_equal(loglik(t1,c(100,1),Model="expo",Rcpp=F),loglik(t1,c(100,1),Model="expo",Rcpp=T))
})


test_that("To ensure the correct output", {
  library(ape)
  t1=rcoal(20)
  x=att(t1)
  expect_true(class(loglik(t1,c(100,1),Model="expo",Rcpp=T))=="numeric")
})
