context("loglikelihood of two methods")

test_that("same loglikelihood", {
  library(ape)
  t1=rcoal(20)
  expect_equal(loglik(t1,c(100,1),Model="expo",Rcpp=F),loglik(t1,c(100,1),Model="expo",Rcpp=T))
})
