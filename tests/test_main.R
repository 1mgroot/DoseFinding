library(testthat)

source("../src/core/main.R")

cat("\n\n==============\nRunning tests for main.R\n==============\n")

test_that("main.R runs without error", {
  expect_true(exists("final_od"))
})

test_that("Final OD is a single value", {
  expect_length(final_od, 1)
})

test_that("Final OD is a dose level", {
  expect_true(final_od %in% trial_config$dose_levels | is.na(final_od))
})

cat("\n==============\nTests for main.R passed\n==============\n")
