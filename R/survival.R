# ---------------------------------------------------------------------------
# Downstream: survival.
#
# Kept deliberately thin. These are convenience wrappers around `survival`
# for the one specific question PGViS is built to ask -- "does the score
# stratify outcome?" -- not a general survival toolkit. Project-specific
# modelling belongs in your analysis scripts, reading pgvis_scores.rds.
# ---------------------------------------------------------------------------

#' Bin a continuous score for survival stratification
#'
#' @param x Numeric score.
#' @param method `"median"`, `"tertile"`, `"quartile"`, or `"cutpoint"`.
#' @param cutpoint Threshold when `method = "cutpoint"`.
#' @return An ordered factor.
#' @export
#' @examples
#' table(pgvis_stratify(rnorm(30), "tertile"))
pgvis_stratify <- function(x, method = c("median", "tertile", "quartile", "cutpoint"),
                           cutpoint = NULL) {
  method <- match.arg(method)
  if (method == "cutpoint") {
    if (is.null(cutpoint)) stop("`cutpoint` is required for method = 'cutpoint'.",
                                call. = FALSE)
    return(factor(ifelse(x > cutpoint, "high", "low"),
                  levels = c("low", "high"), ordered = TRUE))
  }
  probs <- switch(method,
                  median   = c(0, 0.5, 1),
                  tertile  = c(0, 1/3, 2/3, 1),
                  quartile = c(0, 0.25, 0.5, 0.75, 1))
  labs <- switch(method,
                 median   = c("low", "high"),
                 tertile  = c("low", "mid", "high"),
                 quartile = c("Q1", "Q2", "Q3", "Q4"))
  br <- stats::quantile(x, probs, na.rm = TRUE)
  if (anyDuplicated(br)) {
    stop("The score has too many tied values to split into ", length(labs),
         " groups. Use method = 'median' or a manual cutpoint.", call. = FALSE)
  }
  cut(x, breaks = br, labels = labs, include.lowest = TRUE, ordered_result = TRUE)
}

#' Test whether PGViS stratifies survival
#'
#' Runs a log-rank test across score groups and a Cox model on the score,
#' both as a continuous predictor and as the grouped factor. Covariates are
#' supported so you can check the score is not just reporting age or stage.
#'
#' @param scored Output of [pgvis_score()], or any table with the score.
#' @param clinical Table with `sample`, a time column and an event column.
#' @param time_col,event_col Column names. `event_col` must be 0/1 or
#'   logical, with 1/`TRUE` meaning the event occurred.
#' @param score_col Score column.
#' @param covariates Character vector of extra Cox covariates.
#' @param method Stratification for the log-rank test; see [pgvis_stratify()].
#' @return A list with `data`, `fit` (survfit), `logrank`, `cox_continuous`,
#'   `cox_grouped` and a tidy `summary` row.
#' @export
pgvis_survival <- function(scored, clinical,
                           time_col = "time", event_col = "event",
                           score_col = "pgvis",
                           covariates = character(),
                           method = "median") {
  .need("survival", "survival analysis")
  if (!"sample" %in% names(scored))  stop("`scored` needs a `sample` column.", call. = FALSE)
  if (!"sample" %in% names(clinical)) stop("`clinical` needs a `sample` column.", call. = FALSE)
  for (cl in c(time_col, event_col)) {
    if (!cl %in% names(clinical)) stop("`clinical` has no column '", cl, "'.", call. = FALSE)
  }

  df <- dplyr::inner_join(as.data.frame(scored), as.data.frame(clinical),
                          by = "sample")
  n_lost <- nrow(scored) - nrow(df)
  if (n_lost > 0) .msg(n_lost, " scored sample(s) had no clinical record.")

  df$.time  <- as.numeric(df[[time_col]])
  df$.event <- as.integer(df[[event_col]])
  df$.score <- as.numeric(df[[score_col]])
  df <- df[stats::complete.cases(df[, c(".time", ".event", ".score")]), , drop = FALSE]

  if (!nrow(df)) stop("No samples have both a score and complete follow-up.", call. = FALSE)
  if (!all(df$.event %in% c(0L, 1L))) {
    stop("`", event_col, "` must be coded 0/1 (or logical). Found: ",
         paste(utils::head(sort(unique(df$.event))), collapse = ", "), call. = FALSE)
  }
  n_ev <- sum(df$.event)
  if (n_ev < 10L) {
    .msg("Only ", n_ev, " events; Cox estimates will be unstable. ",
         "A rule of thumb is 10 events per covariate.")
  }

  df$.group <- pgvis_stratify(df$.score, method)
  surv <- survival::Surv(df$.time, df$.event)

  fit     <- survival::survfit(surv ~ .group, data = df)
  logrank <- survival::survdiff(surv ~ .group, data = df)
  lr_p    <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1,
                           lower.tail = FALSE)

  rhs <- function(term) paste(c(term, covariates), collapse = " + ")
  cox_c <- survival::coxph(stats::as.formula(paste("surv ~", rhs(".score"))), data = df)
  cox_g <- survival::coxph(stats::as.formula(paste("surv ~", rhs(".group"))), data = df)

  sc <- summary(cox_c)$coefficients[".score", , drop = FALSE]
  ci <- summary(cox_c)$conf.int[".score", , drop = FALSE]

  summary_row <- data.frame(
    n = nrow(df), n_events = n_ev,
    hr_per_sd = unname(exp(sc[, "coef"] * stats::sd(df$.score, na.rm = TRUE))),
    hr = unname(ci[, "exp(coef)"]),
    hr_lower = unname(ci[, "lower .95"]),
    hr_upper = unname(ci[, "upper .95"]),
    cox_p = unname(sc[, ncol(sc)]),
    logrank_p = lr_p,
    stringsAsFactors = FALSE
  )

  list(data = df, fit = fit, logrank = logrank, logrank_p = lr_p,
       cox_continuous = cox_c, cox_grouped = cox_g, summary = summary_row)
}

#' Kaplan-Meier plot of score groups
#'
#' Uses `survminer` when available for a risk table, otherwise falls back to
#' base graphics so the function always works.
#'
#' @param x Output of [pgvis_survival()].
#' @param title Plot title.
#' @param xlab,ylab Axis labels.
#' @return A `ggsurvplot` object, or `NULL` for the base-graphics fallback.
#' @export
pgvis_plot_survival <- function(x, title = NULL,
                                xlab = "Time", ylab = "Survival probability") {
  stopifnot(is.list(x), !is.null(x$fit))
  sub <- sprintf("log-rank p = %s",
                 if (x$logrank_p < 0.001) "< 0.001" else sprintf("%.3f", x$logrank_p))

  if (requireNamespace("survminer", quietly = TRUE)) {
    return(survminer::ggsurvplot(
      x$fit, data = x$data, pval = TRUE, risk.table = TRUE,
      conf.int = TRUE, legend.title = "PGViS",
      title = title, subtitle = sub, xlab = xlab, ylab = ylab,
      ggtheme = ggplot2::theme_bw()))
  }
  .msg("survminer not installed; drawing a base-graphics Kaplan-Meier instead.")
  graphics::plot(x$fit, col = seq_along(x$fit$strata), lwd = 2,
                 xlab = xlab, ylab = ylab,
                 main = paste(c(title, sub), collapse = "\n"))
  graphics::legend("bottomleft", legend = levels(x$data$.group),
                   col = seq_along(x$fit$strata), lwd = 2, bty = "n")
  invisible(NULL)
}
