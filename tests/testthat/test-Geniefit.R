context("Cautious input of Geniefit")

test_that("To ensure the valid input", {
  library(ape)
  t1=rcoal(20)
  expect_warning(expect_error(Geniefit(t1,Model="expo",start=c(100,.1,.1),upper=Inf,lower=-Inf)))
  expect_error(Geniefit(t1,Model="exponential",start=c(100,.1),upper=Inf,lower=0))
  expect_equal(class(Geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0)),"list")
})
