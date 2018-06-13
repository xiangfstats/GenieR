context("Vectorised trajectory")

test_that("a vectorised trajectory", {
  sample<-c(100,0)
  trajectory<-function(x) 100
  expect_error(coalgen_iso(sample, trajectory))
})

test_that("To ensure the correct output", {
  sample1<-cbind(c(9,1,2,1),c(0,.008,.03,.1))
  trajectory<-function(x) exp(10*x)
  example_iso<-coalgen_hetero(sample1, trajectory)
  expect_true(class(example_iso)=="list" && length(example_iso)==2)
})
