library(fcrn)
context("Plotting")


test_that("We get a plot from getSortingPosterior", {
  outt <- getSortingPosterior()
  
  expect_is(outt, "ggplot")
})
