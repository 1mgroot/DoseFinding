library(testthat)
if (basename(getwd()) == "tests") setwd("..")
source("src/core/config.R")
source("src/optimization/threshold_calibration.R")

test_that("threshold calibration forwards the actual allocation method", {
  settings<-default_separate_threshold_settings(FALSE)
  settings$adaptive_allocation_method<-"multinomial"
  config<-create_threshold_trial_config(settings,.45,.75,.6)
  expect_identical(config$adaptive_allocation_method,"multinomial")
})
