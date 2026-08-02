# ---------------------------------------------------------------------------
# Reporting.
#
# Group comparison and the density plot share the same statistics. They live
# here as one implementation each, with the AUC computation switchable via an
# argument.
# ---------------------------------------------------------------------------

#' Compare a score between two groups
#'
#' Computes ROC AUC, Cohen's d and a permutation p-value, and refuses to
#' report them when either group is too small.
#'
#' @param df A table with the score and group columns.
#' @param score_col,group_col Column names.
#' @param min_n Minimum n per group. Below this the comparison is marked
#'   invalid and no statistics are returned, because effect sizes estimated
#'   from a handful of samples are not interpretable.
#' @param n_perm Permutation replicates for the p-value.
#' @param auc Compute ROC AUC (needs `pROC`).
#' @param seed Seed for the permutation test.
#' @return A one-row `data.frame` with the group sizes, statistics, a `valid`
#'   flag and a formatted `label` suitable for annotating a plot.
#' @export
pgvis_compare_groups <- function(df, score_col = "pgvis", group_col = "group",
                                 min_n = 50L, n_perm = 1000L,
                                 auc = TRUE, seed = 42L) {
  g <- droplevels(factor(df[[group_col]]))
  lv <- levels(g)
  blank <- data.frame(group1 = NA_character_, group2 = NA_character_,
                      n1 = NA_integer_, n2 = NA_integer_,
                      auc = NA_real_, cohens_d = NA_real_, perm_p = NA_real_,
                      valid = FALSE, label = NA_character_,
                      stringsAsFactors = FALSE)

  if (length(lv) != 2L) {
    blank$label <- "not a two-group comparison"
    return(blank)
  }
  x <- df[[score_col]][g == lv[1]]
  y <- df[[score_col]][g == lv[2]]
  n1 <- sum(!is.na(x)); n2 <- sum(!is.na(y))

  if (n1 < min_n || n2 < min_n) {
    blank$group1 <- lv[1]; blank$group2 <- lv[2]
    blank$n1 <- n1; blank$n2 <- n2
    blank$label <- sprintf("N = %d vs %d\n(below n = %d, suppressed)", n1, n2, min_n)
    return(blank)
  }

  auc_val <- NA_real_
  if (auc && requireNamespace("pROC", quietly = TRUE)) {
    r <- pROC::roc(response = g, predictor = df[[score_col]],
                   levels = lv, direction = "<", quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(r))
    if (auc_val < 0.5) auc_val <- 1 - auc_val
  }

  d_val <- NA_real_
  if (requireNamespace("effsize", quietly = TRUE)) {
    d_val <- abs(effsize::cohen.d(df[[score_col]], g)$estimate)
  } else {
    sp <- sqrt(((n1 - 1) * stats::var(x, na.rm = TRUE) +
                (n2 - 1) * stats::var(y, na.rm = TRUE)) / (n1 + n2 - 2))
    d_val <- abs(mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / sp
  }

  withr_seed <- function(expr) { old <- .Random.seed_safe(); on.exit(.restore_seed(old))
    set.seed(seed); expr }
  perm_p <- withr_seed({
    obs <- abs(mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE))
    pool <- c(x, y)
    perms <- replicate(n_perm, {
      s <- sample(pool)
      abs(mean(s[seq_len(n1)], na.rm = TRUE) -
          mean(s[(n1 + 1):length(pool)], na.rm = TRUE))
    })
    (sum(perms >= obs) + 1) / (n_perm + 1)   # add-one, so p is never exactly 0
  })

  label <- sprintf("N = %d vs %d (%s vs %s)%s\nCohen's d = %.2f\nPerm. p %s",
                   n1, n2, lv[1], lv[2],
                   if (is.na(auc_val)) "" else sprintf("\nROC AUC = %.2f", auc_val),
                   d_val,
                   if (perm_p < 0.001) "< 0.001" else sprintf("= %.3f", perm_p))

  data.frame(group1 = lv[1], group2 = lv[2], n1 = n1, n2 = n2,
             auc = auc_val, cohens_d = as.numeric(d_val), perm_p = perm_p,
             valid = TRUE, label = label, stringsAsFactors = FALSE)
}

.Random.seed_safe <- function() {
  if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
}
.restore_seed <- function(old) {
  if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
}

#' Density plot of PGViS by group, faceted and annotated
#'
#' @param scored Output of [pgvis_score()], with a `group` column identifying
#'   the comparison arm.
#' @param group_col Column defining the density curves.
#' @param facet_col Optional column to facet by.
#' @param stratum Optional single stratum (e.g. `"EUR"`) to restrict to.
#' @param strata_col Column holding the stratum.
#' @param min_n Minimum group size; panels below it are annotated but not drawn.
#' @param palette Named fill colours.
#' @return A `ggplot` object.
#' @export
pgvis_plot_density <- function(scored, group_col = "group", facet_col = NULL,
                               stratum = NULL, strata_col = "superpopulation",
                               min_n = 50L, palette = NULL) {
  .need("ggplot2", "plotting")
  df <- dplyr::as_tibble(scored)
  if (!is.null(stratum)) df <- dplyr::filter(df, .data[[strata_col]] == stratum)
  df <- dplyr::filter(df, !is.na(.data[[group_col]]), !is.na(.data$pgvis))
  if (!nrow(df)) stop("Nothing left to plot after filtering.", call. = FALSE)

  keys <- if (is.null(facet_col)) {
    stats::setNames(list(df), "all")
  } else {
    split(df, df[[facet_col]])
  }
  stats <- do.call(rbind, lapply(names(keys), function(nm) {
    s <- pgvis_compare_groups(keys[[nm]], "pgvis", group_col, min_n = min_n)
    if (!is.null(facet_col)) s[[facet_col]] <- nm
    s
  }))

  drawable <- df
  if (!is.null(facet_col)) {
    ok <- stats[[facet_col]][stats$valid]
    drawable <- dplyr::filter(df, .data[[facet_col]] %in% ok)
  } else if (!isTRUE(stats$valid)) {
    drawable <- df[0, ]
  }

  stats$x <- max(df$pgvis, na.rm = TRUE)
  stats$y <- Inf

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = drawable,
      ggplot2::aes(x = .data$pgvis, colour = .data[[group_col]],
                   fill = .data[[group_col]]),
      alpha = 0.45, linewidth = 1.1) +
    ggplot2::geom_text(
      data = stats,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 1, vjust = 1.1, size = 3.4, lineheight = 0.95,
      inherit.aes = FALSE) +
    ggplot2::labs(x = "PGViS", y = NULL,
                  colour = NULL, fill = NULL,
                  title = if (is.null(stratum)) NULL else paste("Ancestry:", stratum)) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette) +
      ggplot2::scale_colour_manual(values = palette)
  }
  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)),
                                 scales = "free_y")
  }
  p
}

#' Correlation heatmap of the component metrics
#'
#' Worth looking at before trusting the composite: if two metrics are almost
#' perfectly correlated, PC1 weighting will effectively double-count them.
#'
#' @param df A table containing the metric columns.
#' @param metrics Columns to correlate.
#' @param labels Optional display names.
#' @param title Plot title.
#' @return The correlation matrix, invisibly; called for the plot.
#' @export
pgvis_plot_metric_cor <- function(df,
                                  metrics = c("delta_wi_int", "log_odds_ratio",
                                              "delta_af_alt"),
                                  labels = c(delta_wi_int = "\u0394Wi",
                                             log_odds_ratio = "LOR",
                                             delta_af_alt = "\u0394AF"),
                                  title = "") {
  .need("corrplot", "correlation plots")
  metrics <- intersect(metrics, names(df))
  m <- stats::cor(df[, metrics, drop = FALSE], use = "pairwise.complete.obs")
  nm <- ifelse(metrics %in% names(labels), labels[metrics], metrics)
  dimnames(m) <- list(nm, nm)
  corrplot::corrplot(m, method = "color",
                     col = grDevices::colorRampPalette(
                       c("#2166AC", "white", "#B2182B"))(200),
                     addCoef.col = "black", number.cex = 0.9,
                     tl.col = "black", tl.srt = 45,
                     title = title, mar = c(0, 0, 1, 0))
  invisible(m)
}
