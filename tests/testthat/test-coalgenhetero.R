context("Vectorised trajectory and a n by 2 matrix sample")

test_that("a vectorised trajectory and a n by 2 matrix sample", {
  sample1<-cbind(c(9,1,2,1),c(0,.008,.03,.1))
  trajectory<-function(x) 100
  expect_error(coalgen_hetero(sample1, trajectory))
  expect_error(coalgen_hetero(sample1[,1], Vectorize(trajectory)))
})
