library(testthat)

if (basename(getwd()) == "tests") {
  setwd("..")
}

run_example_smoke_test <- function(filename, expected_files = character()) {
  output_dir <- tempfile(paste0("dosefinding-example-", tools::file_path_sans_ext(filename), "-"))
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  command_output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", file.path("examples", filename)),
    stdout = TRUE,
    stderr = TRUE,
    env = paste0("DOSEFINDING_EXAMPLE_OUTPUT_DIR=", output_dir)
  )

  exit_status <- attr(command_output, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }
  expect_equal(
    exit_status,
    0L,
    info = paste(filename, paste(command_output, collapse = "\n"), sep = "\n")
  )
  if (length(expected_files) > 0) {
    expect_true(
      all(file.exists(file.path(output_dir, expected_files))),
      info = paste("Missing expected output from", filename)
    )
  }
}

test_that("all retained examples run successfully", {
  expected_outputs <- list(
    simple_usage_example.R = c("allocation_plot.png", "posterior_immune_response.png"),
    flat_scenario_demo.R = "flat_scenario_comparison.png",
    plotting_demo.R = c("demo_multi_scenarios.png", "demo_obd_selection.png"),
    logging_control_example.R = character(),
    poc_calibration_demo.R = character()
  )
  retained_examples <- basename(list.files("examples", pattern = "[.]R$"))

  expect_setequal(names(expected_outputs), retained_examples)
  for (filename in names(expected_outputs)) {
    run_example_smoke_test(filename, expected_outputs[[filename]])
  }
})
