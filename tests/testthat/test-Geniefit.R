context("Cautious input of Geniefit")

test_that("To ensure the valid input", {
  library(ape)
  t1=rcoal(20)
  expect_error(Geniefit(t1,Model="expo",start=c(100,.1,.1),upper=Inf,lower=-Inf))
  expect_error(Geniefit(t1,Model="exponential",start=c(100,.1),upper=Inf,lower=0))
  expect_equal(class(Geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0)),"list")
})

test_that("numbers of parameters", {
  library(ape)
  t1=rcoal(20)
  expect_error(Geniefit(t1,Model="const",start=c(100,.1),upper=Inf,lower=0))
  expect_error(Geniefit(t1,Model="expo",start=c(100),upper=Inf,lower=0))
  expect_error(Geniefit(t1,Model="log",start=c(100,.1),upper=Inf,lower=0))
  expect_error(Geniefit(t1,Model="step",start=c(100,.1),upper=Inf,lower=0))
  expect_error(Geniefit(t1,Model="plog",start=c(100,.1),upper=Inf,lower=0))
  expect_error(Geniefit(t1,Model="pexpan",start=c(100,.1),upper=Inf,lower=0))
  expect_error(Geniefit(t1,Model="expan",start=c(100,.1),upper=Inf,lower=0))
})

test_that("correct outputs", {
  library(ape)
  t1=rcoal(20)
  expect_true(class(Geniefit(t1,Model="const",start=c(100),upper=Inf,lower=0))=="list"  )
  expect_true(class(Geniefit(t1,Model="expo",start=c(100,.1),upper=Inf,lower=0))=="list"  )
  expect_true(class(Geniefit(t1,Model="log",start=c(100,.1,.1),upper=Inf,lower=0))=="list"  )
  expect_true(class(Geniefit(t1,Model="step",start=c(100,.1,.1),upper=Inf,lower=0))=="list"  )
  expect_true(class(Geniefit(t1,Model="plog",start=c(100,.1,.1),upper=Inf,lower=0))=="list"  )
  expect_true(class(Geniefit(t1,Model="pexpan",start=c(100,.1,.1),upper=Inf,lower=0))=="list"  )
  expect_true(class(Geniefit(t1,Model="expan",start=c(100,.1,.1),upper=Inf,lower=0))=="list"   )
})
