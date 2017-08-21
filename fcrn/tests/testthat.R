# see https://github.com/hadley/testthat/issues/144
Sys.setenv(R_TESTS = '')

library(testthat)
library(fcrn)

test_check('fcrn')