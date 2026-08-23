# tests/testthat/test-survival.R

library(testthat)

test_that("km_fit_and_plot returns plot and surv objects", {
  # simulate small dataset
  set.seed(1)
  n <- 50
  df <- data.frame(OS_time = rexp(n, 0.1), OS_event = rbinom(n,1,0.6), Biomarker_Status = sample(c('Low','High'), n, replace = TRUE))
  res <- km_fit_and_plot(df, time = 'OS_time', event = 'OS_event', group = 'Biomarker_Status')
  expect_true(!is.null(res$plot))
  expect_true(inherits(res$survfit, 'survfit'))
})
