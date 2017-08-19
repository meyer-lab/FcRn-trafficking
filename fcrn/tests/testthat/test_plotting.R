library(fcrn)
context("Plotting")


test_that("We get a plot from getSortingPosterior", {
  outt <- getSortingPosterior()
  
  expect_is(outt, "ggplot")
})

test_that("We get a plot from plot_halfls", {
  outt <- plot_halfls("diff")
  
  expect_is(outt, "ggplot")
})
