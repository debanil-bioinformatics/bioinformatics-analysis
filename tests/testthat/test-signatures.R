# tests/testthat/test-signatures.R

library(testthat)

test_that("score_by_mean returns expected dimensions", {
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0('Gene', 1:10)
  scores <- score_by_mean(mat, c('Gene1','Gene2','Gene3'))
  expect_equal(length(scores), ncol(mat))
})
