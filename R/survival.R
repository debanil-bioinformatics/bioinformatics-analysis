# survival.R - reusable survival analysis helpers

# km_fit_and_plot: fit Kaplan-Meier and optionally Cox model, return plot and objects
km_fit_and_plot <- function(data, time, event, group, group_levels = NULL,
                            max_followup = NULL, conf.int = FALSE,
                            annotate = TRUE, title = NULL, xlab = "Time",
                            ylab = "Survival Probability", ...) {
  stopifnot(all(c(time, event, group) %in% colnames(data)))
  d <- data
  if (!is.null(group_levels)) d[[group]] <- factor(d[[group]], levels = group_levels)
  if (!is.null(max_followup)) d$surv_time <- pmin(d[[time]], max_followup) else d$surv_time <- d[[time]]
  d$event_status <- d[[event]]
  surv_obj <- survival::Surv(time = d$surv_time, event = d$event_status)
  fit <- survival::survfit(surv_obj ~ d[[group]], data = d)
  cox <- tryCatch(survival::coxph(surv_obj ~ d[[group]], data = d), error = function(e) NULL)

  p <- survminer::ggsurvplot(
    fit, data = d, conf.int = conf.int, risk.table = FALSE,
    title = title, xlab = xlab, ylab = ylab, legend.title = group,
    legend.labs = if (!is.null(group_levels)) group_levels else levels(d[[group]]), ...
  )

  if (annotate && !is.null(cox)) {
    s <- summary(cox)
    hr <- round(s$coefficients[1, "exp(coef)"], 2)
    pv <- format.pval(s$coefficients[1, "Pr(>|z|)"], digits = 3)
    ann <- paste0("HR=", hr, ", p=", pv)
    p$plot <- p$plot + ggplot2::annotate("text", x = Inf, y = Inf, label = ann,
                                       hjust = 1.1, vjust = 1.5, size = 3)
  }

  list(plot = p$plot, survfit = fit, cox = cox)
}

# cox_multivariable: fit multivariable Cox model given covariates (character vector)
cox_multivariable <- function(data, time, event, covariates) {
  stopifnot(all(c(time, event) %in% colnames(data)))
  vars <- unique(covariates)
  stopifnot(all(vars %in% colnames(data)))
  form <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", paste(vars, collapse = " + ")))
  model <- survival::coxph(form, data = data)
  summary(model)
}
