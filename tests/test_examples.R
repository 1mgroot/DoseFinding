library(testthat)

if (basename(getwd()) == "tests") {
  setwd("..")
}

test_that("simple usage example runs successfully", {
  output_dir <- tempfile("dosefinding-simple-example-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  command_output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "examples/simple_usage_example.R"),
    stdout = TRUE,
    stderr = TRUE,
    env = paste0("DOSEFINDING_EXAMPLE_OUTPUT_DIR=", output_dir)
  )

  exit_status <- attr(command_output, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }
  expect_equal(exit_status, 0L, info = paste(command_output, collapse = "\n"))
  expect_true(file.exists(file.path(output_dir, "allocation_plot.png")))
  expect_true(file.exists(file.path(output_dir, "posterior_immune_response.png")))
})
