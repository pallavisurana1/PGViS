# ---------------------------------------------------------------------------
# Sample-level scoring.
#
# Variant-level metrics are aggregated to one value per sample, standardised
# against an ancestry-matched reference cohort, and combined into a single
# composite score. Ancestry matching is not cosmetic: allele frequency
# divergence from a reference panel is strongly confounded by how genetically
# distant the sample is from that panel, so the standardisation is done
# *within* superpopulation.
# ---------------------------------------------------------------------------

#' Expand a carrier-sample column into one row per sample
#'
#' Cohort BEDs store carriers as a delimited string. Scoring needs one row per
#' (sample, variant).
#'
#' @param df A variant table.
#' @param samples_col Column holding the delimited sample IDs.
#' @param sample_col Name for the new per-sample column.
#' @param sep Regex separating sample IDs.
#' @return A long table with `n_samples` recording the original carrier count.
#' @export
#' @examples
#' pgvis_explode_carriers(
#'   data.frame(var_site = "chr1:1-2(A>G)", carrier_samples = "S1;S2,S3")
#' )
pgvis_explode_carriers <- function(df,
                                   samples_col = "carrier_samples",
                                   sample_col = "sample",
                                   sep = "\\s*[;,]\\s*") {
  if (!samples_col %in% names(df)) {
    stop("No column named '", samples_col, "'.", call. = FALSE)
  }
  toks <- stringr::str_split(as.character(df[[samples_col]]), sep)
  toks <- lapply(toks, function(v) {
    v <- trimws(v)
    v[nzchar(v) & v != "." & !is.na(v)]
  })
  df$n_samples <- lengths(toks)
  df[[sample_col]] <- toks
  out <- tidyr::unnest(dplyr::as_tibble(df), cols = dplyr::all_of(sample_col))
  dplyr::filter(out, !is.na(.data[[sample_col]]), .data[[sample_col]] != "")
}

#' Aggregate variant metrics to one row per sample
#'
#' @param df A long (sample x variant) table.
#' @param metrics Metric columns to average.
#' @param by Grouping columns; `sample` plus whatever identifies the cohort
#'   and ancestry.
#' @param fun Aggregation function.
#' @return One row per group, with `n_variants` and the aggregated metrics.
#' @export
pgvis_aggregate_samples <- function(df,
                                    metrics = c("delta_wi_int", "log_odds_ratio",
                                                "delta_af_alt"),
                                    by = c("sample", "superpopulation", "cohort"),
                                    fun = mean) {
  by <- intersect(by, names(df))
  metrics <- intersect(metrics, names(df))
  if (!length(metrics)) stop("None of the requested metric columns are present.",
                             call. = FALSE)
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_variants = dplyr::n(),
      dplyr::across(dplyr::all_of(metrics), \(x) fun(x, na.rm = TRUE)),
      .groups = "drop"
    )
}

#' Compute the PGViS composite score
#'
#' Three steps:
#'
#' 1. **Robust standardisation.** Each metric is centred on the median and
#'    scaled by the median absolute deviation *of the reference cohort within
#'    the same superpopulation*. Median/MAD rather than mean/SD because the
#'    metric distributions have heavy right tails.
#' 2. **Weighting.** Weights are the absolute PC1 loadings of the scaled
#'    z-scores, normalised to sum to one. PC1 is the axis of shared variation
#'    across the three metrics, so this weights each metric by how much it
#'    contributes to the common signal rather than by a hand-picked constant.
#' 3. **Combination.** The score is the weighted sum of the z-scores.
#'
#' @param df One row per sample, from [pgvis_aggregate_samples()].
#' @param metrics Metric columns to combine.
#' @param reference_group Value of `group_col` identifying the reference
#'   cohort (e.g. `"1000Genomes"`).
#' @param group_col Column holding the cohort label.
#' @param strata Column defining the ancestry strata. Set to `NULL` to
#'   standardise against the whole reference cohort, which is only appropriate
#'   for an ancestry-homogeneous dataset.
#' @param weights Optional named numeric vector to override the PCA-derived
#'   weights. Useful for scoring a new cohort with weights frozen from a
#'   previous run; see [pgvis_weights()].
#' @return A table with one row per sample: the per-metric z-scores, the
#'   composite `pgvis` score, and the weights as attribute `"weights"`.
#' @export
pgvis_score <- function(df,
                        metrics = c("delta_wi_int", "log_odds_ratio", "delta_af_alt"),
                        reference_group = "1000Genomes",
                        group_col = "cohort",
                        strata = "superpopulation",
                        weights = NULL) {

  metrics <- intersect(metrics, names(df))
  if (length(metrics) < 2L) {
    stop("At least two metric columns are needed to compute a composite score.",
         call. = FALSE)
  }
  if (!group_col %in% names(df)) stop("No column '", group_col, "'.", call. = FALSE)
  if (!reference_group %in% df[[group_col]]) {
    stop("Reference group '", reference_group, "' not present in `", group_col,
         "`. The score is defined relative to a reference cohort; without one there is nothing to standardise against.",
         call. = FALSE)
  }
  if (!is.null(strata) && !strata %in% names(df)) {
    stop("No stratum column '", strata, "'. Pass strata = NULL to disable ancestry matching.",
         call. = FALSE)
  }

  ref <- dplyr::filter(df, .data[[group_col]] == reference_group)
  if (!is.null(strata)) ref <- dplyr::group_by(ref, .data[[strata]])
  ref <- dplyr::summarise(
    ref,
    dplyr::across(dplyr::all_of(metrics),
                  list(mu = \(x) stats::median(x, na.rm = TRUE),
                       sg = \(x) stats::mad(x, na.rm = TRUE)),
                  .names = "{.col}__{.fn}"),
    .groups = "drop"
  )

  # A zero MAD makes the z-score undefined; fall back to the SD.
  for (m in metrics) {
    sg <- paste0(m, "__sg")
    bad <- !is.finite(ref[[sg]]) | ref[[sg]] == 0
    if (any(bad)) {
      .msg("MAD is zero for '", m, "' in ", sum(bad),
           " stratum/strata; falling back to SD there.")
      fb <- stats::sd(df[[m]], na.rm = TRUE)
      ref[[sg]][bad] <- if (is.finite(fb) && fb > 0) fb else 1
    }
  }

  out <- if (is.null(strata)) {
    dplyr::bind_cols(df, ref[rep(1L, nrow(df)), , drop = FALSE])
  } else {
    dplyr::left_join(df, ref, by = strata)
  }

  zcols <- character(0)
  for (m in metrics) {
    z <- paste0("z_", m)
    out[[z]] <- (out[[m]] - out[[paste0(m, "__mu")]]) / out[[paste0(m, "__sg")]]
    zcols <- c(zcols, z)
  }

  scaled <- as.data.frame(lapply(out[zcols], \(x) as.numeric(scale(x))))
  names(scaled) <- zcols

  if (is.null(weights)) {
    weights <- pgvis_weights(scaled)
  } else {
    miss <- setdiff(zcols, names(weights))
    if (length(miss)) stop("`weights` is missing: ", paste(miss, collapse = ", "),
                           call. = FALSE)
    weights <- weights[zcols] / sum(weights[zcols])
  }

  out$pgvis <- as.numeric(as.matrix(scaled[zcols]) %*% weights[zcols])

  keep <- c(intersect(c("sample", strata, group_col, "n_variants"), names(out)),
            metrics, zcols, "pgvis")
  res <- out[, keep, drop = FALSE]
  attr(res, "weights") <- weights
  attr(res, "reference_group") <- reference_group
  res
}

#' Derive metric weights from the first principal component
#'
#' Weights are the absolute PC1 loadings, normalised to sum to one. Exported
#' so that weights fitted on a discovery cohort can be frozen and reused when
#' scoring a validation cohort, rather than being silently refitted.
#'
#' @param scaled A data frame of scaled z-scores, one column per metric.
#' @return A named numeric vector summing to 1.
#' @export
pgvis_weights <- function(scaled) {
  cc <- stats::complete.cases(scaled)
  if (sum(cc) < 3L) {
    .msg("Too few complete cases for PCA; falling back to equal weights.")
    w <- rep(1 / ncol(scaled), ncol(scaled))
    names(w) <- names(scaled)
    return(w)
  }
  fit <- stats::prcomp(scaled[cc, , drop = FALSE], center = TRUE, scale. = TRUE)
  load <- abs(fit$rotation[, "PC1"])
  w <- load / sum(load)
  .msg("PC1-derived weights: ",
       paste(sprintf("%s=%.3f", names(w), w), collapse = "  "),
       "  (PC1 explains ",
       sprintf("%.1f%%", 100 * fit$sdev[1]^2 / sum(fit$sdev^2)), " of variance)")
  w
}

#' Map self-reported race to a 1000 Genomes superpopulation
#'
#' A convenience mapping for TCGA/CPTAC-style clinical fields. Self-reported
#' race is a social variable and only a rough proxy for genetic ancestry;
#' where genotype data are available, prefer ancestry inferred from principal
#' components or a tool such as `peddy`. Unmapped values become `NA` and those
#' samples are excluded from scoring rather than pooled into a default group.
#'
#' @param race Character vector of self-reported race.
#' @param mapping Named character vector; names are lower-cased race strings.
#' @return A factor with levels `AFR`, `AMR`, `EAS`, `EUR`, `SAS`.
#' @export
#' @examples
#' pgvis_map_superpopulation(c("white", "asian", "not reported"))
pgvis_map_superpopulation <- function(race, mapping = NULL) {
  default <- c(
    "black or african american"        = "AFR",
    "african american"                 = "AFR",
    "american indian or alaska native" = "AMR",
    "asian"                            = "EAS",
    "native hawaiian or other pacific islander" = "EAS",
    "white"                            = "EUR"
  )
  mapping <- mapping %||% default
  key <- tolower(trimws(as.character(race)))
  out <- unname(mapping[key])
  factor(out, levels = c("AFR", "AMR", "EAS", "EUR", "SAS"))
}
