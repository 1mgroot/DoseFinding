library(testthat)
if (basename(getwd()) == "tests") setwd("..")
source("src/core/main.R")

test_that("missing strata receive actual Beta(1,1) draws and variance", {
  set.seed(42)
  x <- data.frame(d=c(1,1,2),Y_I=c(0,1,0),r=c(0,1,0),n=c(2,1,2))
  p <- add_beta_variance(simulate_beta_posterior(x,n_sims=8))
  expect_warning(result <- apply_biviso_on_matrix(p), NA)
  missing <- result[result$d==2 & result$Y_I==1,]
  expect_equal(missing$var_post, 1/12)
  expect_equal(missing$alpha_post,1)
  expect_equal(missing$beta_post,1)
  expect_length(missing$samples[[1]],8)
  expect_true(all(is.finite(missing$samples[[1]])))
  expect_equal(result$samples[1:3],p$samples)
})

test_that("entire absent immune group and single dose are supported", {
  x <- data.frame(d=1:2,Y_I=0,r=0,n=3)
  p <- add_beta_variance(simulate_beta_posterior(x,n_sims=8))
  expect_warning(result <- apply_biviso_on_matrix(p), NA)
  expect_equal(nrow(result),4)
  expect_true(all(lengths(result$samples)==8))
  one <- apply_biviso_on_matrix(p[1,,drop=FALSE])
  expect_equal(nrow(one),2)
  expect_true(all(one$samples_pava[[1]]<=one$samples_pava[[2]]))
})

test_that("marginalization aligns by dose labels rather than row positions", {
  x <- data.frame(d=c(10,10,20,20),Y_I=rep(0:1,2),
    samples_pava=I(list(rep(.1,3),rep(.8,3),rep(.2,3),rep(.9,3))))
  immune <- data.frame(d=c(20,10),samples_pava=I(list(rep(.8,3),rep(.1,3))))
  result <- compute_marginal_probability(x,immune)
  expect_equal(result$d,c(20,10))
  expect_equal(result$marginal_prob,c(.76,.17))
  expect_equal(result,compute_marginal_probability(x[c(4,2,1,3),],immune))
  x$d[x$d==20] <- 30
  expect_error(compute_marginal_probability(x,immune),"Dose labels must match")
})

test_that("malformed draw counts cannot be silently recycled", {
  x <- data.frame(d=rep(1:2,each=2),Y_I=rep(0:1,2),r=0,n=2)
  p <- add_beta_variance(simulate_beta_posterior(x,n_sims=8))
  p$samples[[1]] <- p$samples[[1]][1:4]
  expect_error(apply_biviso_on_matrix(p),"equal draw counts")
})

test_that("single-dose immune adjustment preserves the draw dimension", {
  p <- add_beta_variance(simulate_beta_posterior(data.frame(d=10,r=1,n=3),n_sims=8))
  adjusted <- apply_pava_on_samples(p)
  expect_equal(nrow(adjusted),1)
  expect_identical(adjusted$samples_pava[[1]],p$samples[[1]])
})
