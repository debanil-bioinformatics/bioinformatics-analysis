# regression.R - regression helpers

# tidy_cox_summary: return data.frame of HR, CI, p for coxph model
tidy_cox_summary <- function(cox_model) {
  stopifnot(inherits(cox_model, "coxph"))
  s <- summary(cox_model)
  coefs <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  res <- data.frame(
    term = rownames(coefs),
    HR = ci$`exp(coef)`,
    HR_low = ci$`lower .95`,
    HR_high = ci$`upper .95`,
    p = coefs$`Pr(>|z|)`
  )
  res
}
