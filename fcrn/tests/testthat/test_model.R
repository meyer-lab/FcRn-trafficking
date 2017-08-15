library(fcrn)
context("Model construction")


test_that("rstan can compile the model", {
  file_path <- system.file("extdata", "model.stan", package = "fcrn")
  
  output <- rstan::stan_model(file_path, auto_write = T, save_dso = T, verbose = F)
  
  expect_is(output, "stanmodel")
})

test_that("Jacobian is build successfully", {
  p <- list(Q = 1, Vin = 1, sortF = 0.9, releaseF = 0.8, Vp = 1, Qu = 1)
  
  output <- jacFunc(p)
  
  expect_is(output, "matrix")
  
  
  outt <- halfl_fcrn(p)
  
  expect_lt(outt, 1.2)
  expect_gt(outt, 1.0)
})