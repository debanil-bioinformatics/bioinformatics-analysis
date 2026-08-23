# tests/testthat/test-correlation.R

library(testthat)

test_that("gene_correlation computes expected Pearson", {
  x <- 1:10
  y <- 2*(1:10) + rnorm(10,0,1)
  res <- gene_correlation(x,y, method = 'pearson')
  expect_true(abs(res$r) > 0.8)
  expect_true(res$p.value < 0.05)
})
