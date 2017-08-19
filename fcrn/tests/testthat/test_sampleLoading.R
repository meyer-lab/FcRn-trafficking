library(fcrn)
context("Sample loading")


test_that("we can load the different model samples", {
  outt <- loadsample("diff")
  
  expect_is(outt, "stanfit")
  
  outt <- loadsample("scarlette")
  
  expect_is(outt, "stanfit")
  
  outt <- loadsample("marlene")
  
  expect_is(outt, "stanfit")
  
  outt <- loadsample("diff", as.df = T)
  
  expect_is(outt, "data.frame")
  
  outt <- loadsample("scarlette", as.df = T)
  
  expect_is(outt, "data.frame")
  
  outt <- loadsample("marlene", as.df = T)
  
  expect_is(outt, "data.frame")
})


test_that("we can load all at once", {
  outt <- loadAll()
  
  expect_is(outt, "data.frame")
})