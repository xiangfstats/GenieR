context("Vectorised trajectory")

test_that("a vectorised trajectory", {
  sample<-c(100,0)
  trajectory<-function(x) 100
  expect_error(coalgen_iso(sample, trajectory))
})
