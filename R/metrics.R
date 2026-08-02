# ---------------------------------------------------------------------------
# Per-variant metrics.
#
# Three quantities are computed for every (region, variant) pair:
#
#   log_odds_ratio  - how much the model's confidence in the regulatory
#                     element changes, on a log2-odds scale
#   score_change    - a magnitude-weighted probability difference
#   delta_wi        - divergence between the cohort allele frequencies and
#                     the matched population frequencies
#
# `log_odds_ratio` and `score_change` were originally defined twice in R and
# once in Python, with `log_odds_ratio` and `log_odd_ratio` used
# interchangeably as column names. One definition each, here.
# ---------------------------------------------------------------------------

#' Log-odds ratio between reference and alternate model probabilities
#'
#' \deqn{\mathrm{LOR} = \log_2\frac{p_{ref}}{1-p_{ref}} - \log_2\frac{p_{alt}}{1-p_{alt}}}
#'
#' Positive values mean the variant reduces the model's confidence that the
#' sequence is a functional element, i.e. the variant is disruptive.
#'
#' @param p_ref,p_alt Model probabilities on the reference and alternate
#'   sequences, in (0, 1).
#' @param eps Probabilities are clamped to `[eps, 1 - eps]` so that saturated
#'   predictions give a large finite value instead of `Inf`.
#' @return Numeric vector.
#' @export
#' @examples
#' log_odds_ratio(0.98, 0.02)
log_odds_ratio <- function(p_ref, p_alt, eps = 1e-6) {
  p_ref <- pmin(pmax(p_ref, eps), 1 - eps)
  p_alt <- pmin(pmax(p_alt, eps), 1 - eps)
  log2(p_ref / (1 - p_ref)) - log2(p_alt / (1 - p_alt))
}

#' Magnitude-weighted probability change
#'
#' \deqn{\Delta = (p_{alt} - p_{ref}) \cdot \max(p_{ref}, p_{alt})}
#'
#' Down-weights changes between two low probabilities, where the model is not
#' confident either way.
#'
#' @param p_ref,p_alt Model probabilities.
#' @return Numeric vector; negative for disruptive variants.
#' @export
#' @examples
#' score_change(0.9, 0.1)
score_change <- function(p_ref, p_alt) {
  (p_alt - p_ref) * pmax(p_ref, p_alt)
}

#' Allele-frequency divergence (delta-Wi)
#'
#' Compares the ratio of cohort to population reference-allele frequency
#' against the same ratio for the alternate allele:
#'
#' \deqn{\Delta W_i = \left| \frac{f^{coh}_{ref}}{f^{pop}_{ref}} - \frac{f^{coh}_{alt}}{f^{pop}_{alt}} \right|}
#'
#' A variant that is enriched in the cohort relative to the matched
#' population gives a large value. Because the raw statistic is heavily
#' right-skewed and unbounded, it is rank-based inverse-normal transformed
#' before being combined with the other metrics.
#'
#' @param df A table with cohort and population frequency columns.
#' @param cohort_ref,cohort_alt,pop_ref,pop_alt Column names.
#' @param pseudocount Added to every frequency to avoid division by zero.
#' @param transform Apply the rank-based inverse normal transform. Note this
#'   is computed *within the supplied table*, so pass the full variant set
#'   rather than calling this per group.
#' @param keep_raw Also return the untransformed `delta_wi`.
#' @return `df` with `delta_af_alt`, `delta_af_ref` and `delta_wi_int` added.
#' @export
#' @examples
#' d <- data.frame(af_ref = 0.7, af_alt = 0.3,
#'                 pop_af_ref = 0.9, pop_af_alt = 0.1)
#' pgvis_delta_wi(d, transform = FALSE, keep_raw = TRUE)$delta_wi
pgvis_delta_wi <- function(df,
                           cohort_ref = "af_ref",
                           cohort_alt = "af_alt",
                           pop_ref    = "pop_af_ref",
                           pop_alt    = "pop_af_alt",
                           pseudocount = 1e-6,
                           transform = TRUE,
                           keep_raw = FALSE) {

  if (!pop_ref %in% names(df)) {
    if (!pop_alt %in% names(df)) {
      stop("Need at least one population frequency column.", call. = FALSE)
    }
    df[[pop_ref]] <- 1 - df[[pop_alt]]
  }
  for (cl in c(cohort_ref, cohort_alt, pop_alt)) {
    if (!cl %in% names(df)) stop("Missing column: ", cl, call. = FALSE)
  }

  cr <- df[[cohort_ref]] + pseudocount
  ca <- df[[cohort_alt]] + pseudocount
  pr <- df[[pop_ref]]    + pseudocount
  pa <- df[[pop_alt]]    + pseudocount

  df$delta_af_alt <- df[[cohort_alt]] - df[[pop_alt]]
  df$delta_af_ref <- df[[cohort_ref]] - df[[pop_ref]]

  raw <- abs((cr / pr) - (ca / pa))
  if (keep_raw) df$delta_wi <- raw
  df$delta_wi_int <- if (transform) rank_inverse_normal(raw) else raw
  df
}

#' Rank-based inverse normal transform
#'
#' Maps a vector onto standard normal quantiles via its ranks (Blom-style
#' with an offset of 0.5). `NA`s are preserved and excluded from the ranking.
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length.
#' @export
#' @examples
#' round(rank_inverse_normal(c(1, 5, 2, 100)), 3)
rank_inverse_normal <- function(x) {
  ok <- !is.na(x)
  n <- sum(ok)
  if (n < 2L) return(rep(NA_real_, length(x)))
  out <- rep(NA_real_, length(x))
  out[ok] <- stats::qnorm((rank(x[ok]) - 0.5) / n)
  out
}

#' Attach model predictions and compute all per-variant metrics
#'
#' The one call that takes a variant table plus model probabilities and
#' returns everything the scoring stage needs.
#'
#' @param variants Variant table with cohort and population frequencies.
#' @param predictions Table with `feature`, `func_site`, `var_site`,
#'   `p_ref`, `p_alt`. See [pgvis_read_predictions()].
#' @param cfg A [pgvis_config()]; supplies `pseudocount`.
#' @param by Join keys.
#' @return A table with `log_odds_ratio`, `score_change`, `delta_af_alt`,
#'   `delta_af_ref` and `delta_wi_int`.
#' @export
pgvis_variant_metrics <- function(variants, predictions, cfg = NULL,
                                  by = c("feature", "func_site", "var_site")) {
  miss <- setdiff(by, names(predictions))
  if (length(miss)) stop("Predictions lack join key(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)

  out <- dplyr::inner_join(variants, predictions, by = by)
  n_drop <- nrow(variants) - nrow(out)
  if (n_drop > 0) .msg(sprintf("%d variant(s) had no model prediction and were dropped.", n_drop))
  if (!nrow(out)) {
    stop("No variants matched a prediction. Check that `feature` labels and coordinate keys agree between the two tables.",
         call. = FALSE)
  }

  out$log_odds_ratio <- log_odds_ratio(out$p_ref, out$p_alt)
  out$score_change   <- score_change(out$p_ref, out$p_alt)
  pseudocount <- if (is.null(cfg)) 1e-6 else cfg$thresholds$pseudocount
  pgvis_delta_wi(out, pseudocount = pseudocount)
}

#' Filter to disruptive variants
#'
#' Keeps variants where the model is confident the reference sequence is a
#' functional element and confident that the mutated sequence is not, and
#' where the variant is carried by enough samples to be worth scoring.
#'
#' @param df Output of [pgvis_variant_metrics()].
#' @param cfg A [pgvis_config()] supplying the thresholds.
#' @param ref_prob_min,alt_prob_max,min_carrier_frac Override `cfg`.
#' @param drop_freq_sources Population sources to exclude. Defaults to
#'   `cfg$thresholds$drop_freq_sources`. Excluding `"1000G"` matters when the
#'   comparison cohort is itself 1000 Genomes.
#' @return The filtered table, with an attribute `"pgvis_filter"` recording
#'   how many rows each criterion removed.
#' @export
pgvis_filter_disruptive <- function(df, cfg = NULL,
                                    ref_prob_min = NULL,
                                    alt_prob_max = NULL,
                                    min_carrier_frac = NULL,
                                    drop_freq_sources = NULL) {
  th <- if (is.null(cfg)) list() else cfg$thresholds
  ref_prob_min      <- ref_prob_min      %||% th$ref_prob_min      %||% 0.5
  alt_prob_max      <- alt_prob_max      %||% th$alt_prob_max      %||% 0.5
  min_carrier_frac  <- min_carrier_frac  %||% th$min_carrier_frac  %||% 0
  drop_freq_sources <- drop_freq_sources %||% th$drop_freq_sources %||% character()

  audit <- c(input = nrow(df))

  df <- dplyr::filter(df, .data$p_ref > ref_prob_min, .data$p_alt <= alt_prob_max)
  audit["after_model_prob"] <- nrow(df)

  if (length(drop_freq_sources) && "freq_source" %in% names(df)) {
    df <- dplyr::filter(df, !.data$freq_source %in% drop_freq_sources)
    audit["after_freq_source"] <- nrow(df)
  }
  if (min_carrier_frac > 0 && "n_carrier" %in% names(df) && nrow(df)) {
    cut <- min_carrier_frac * max(df$n_carrier, na.rm = TRUE)
    df <- dplyr::filter(df, .data$n_carrier >= cut)
    audit["after_carrier_frac"] <- nrow(df)
  }

  attr(df, "pgvis_filter") <- audit
  df
}

#' Count unique sites in a variant table
#'
#' @param df A variant table.
#' @param label A label for the returned row.
#' @return A one-row `data.frame`.
#' @export
pgvis_count_sites <- function(df, label = "dataset") {
  pick <- function(nm) if (nm %in% names(df)) length(unique(df[[nm]])) else NA_integer_
  data.frame(
    dataset        = label,
    functional_sites = pick("func_site"),
    variant_sites    = pick("var_site"),
    dbsnp_ids        = pick("rsid"),
    rows             = nrow(df),
    stringsAsFactors = FALSE
  )
}
