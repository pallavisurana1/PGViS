# ---------------------------------------------------------------------------
# Sample-level scoring.
#
# The PGViS score (v0.2.0) is built in four steps:
#
#   1. each variant feature is converted to a percentile within the variant
#      list of its own dataset            -> pgvis_variant_score()
#   2. the percentiles of |delta Wi| and the log odds ratio are combined by
#      their geometric mean to give one score per variant
#                                         -> pgvis_variant_score()
#   3. a sample's score is the root mean square of the scores of the variants
#      it carries                        -> pgvis_sample_scores()
#   4. the sample score is standardised against the reference cohort scored the
#      same way, *within* superpopulation  -> pgvis_score()
#
# In formulas, for a dataset with variant list U and a sample i carrying C_i:
#
#   p_fv    = #{u in U : x_fu <= x_fv} / |U|,   f in {|dWi|, LOR}        (1)
#   s_v     = (p_|dWi|,v * p_LOR,v)^(1/2)                                (2)
#   S_i     = sqrt( (1 / |C_i|) * sum_{v in C_i} s_v^2 )                 (3)
#   PGViS_i = (S_i - mean_ref) / SD_ref                                  (4)
#
# The alternate-allele-frequency shift (delta AF) is deliberately NOT part of
# the score: it is the criterion used to select variants (|delta AF| >= 0.1),
# so tumour delta AF is enriched by construction while the reference panel's
# is near zero, and using it again would make the comparison circular. The
# cumulative sum of variant scores is also avoided because it simply tracks
# how many variants a cohort's variant list contains.
#
# The pre-0.2.0 PCA-weighted composite is still available as
# `pgvis_score(method = "pca")`.
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

#' Score each variant: percentile within its dataset, then geometric mean
#'
#' Implements Eq. 1 and 2 of the PGViS score. Every feature is converted to its
#' percentile among the variants of the *same dataset* (the variant list that
#' passed selection in that cohort), so features with different units and
#' ranges are comparable and the score does not depend on how many variants a
#' cohort contains. The percentiles are then combined by their geometric mean.
#'
#' The default features are the magnitude of the regulatory interaction change
#' (`|delta_wi_int|`) and the DNABERT-derived log odds ratio
#' (`log_odds_ratio`). `delta_af_alt` is intentionally not a default: it is
#' the variant-selection criterion, so including it is circular.
#'
#' @param long A long (sample x variant) table, e.g. from
#'   [pgvis_explode_carriers()]. Reference and test cohorts should be stacked
#'   in one table and distinguished by `universe_col`.
#' @param features Feature columns to rank and combine.
#' @param variant_col Column identifying the variant (default `"var_site"`).
#' @param universe_col Column defining the dataset whose variant list sets the
#'   percentiles (default `"cohort"`). Each distinct value gets its own list.
#' @param abs_features Features ranked by absolute value.
#' @return `long` with one percentile column `p_<feature>` per feature and the
#'   combined `variant_score`. The number of variants in each dataset's list is
#'   attached as attribute `"universe_size"`.
#' @export
#' @examples
#' long <- data.frame(
#'   cohort = "tumour", var_site = c("v1", "v2", "v3"),
#'   delta_wi_int = c(0.1, -0.5, 0.9), log_odds_ratio = c(1, 3, 2)
#' )
#' pgvis_variant_score(long)
pgvis_variant_score <- function(long,
                                features = c("delta_wi_int", "log_odds_ratio"),
                                variant_col = "var_site",
                                universe_col = "cohort",
                                abs_features = "delta_wi_int") {
  miss <- setdiff(c(variant_col, universe_col, features), names(long))
  if (length(miss)) {
    stop("`long` is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!length(features)) stop("Give at least one feature.", call. = FALSE)
  if ("delta_af_alt" %in% features) {
    .msg("delta_af_alt is the variant-selection criterion; including it in the ",
         "score makes the tumour-versus-reference comparison circular.")
  }

  vals <- lapply(features, function(f) {
    x <- suppressWarnings(as.numeric(long[[f]]))
    if (f %in% abs_features) abs(x) else x
  })
  names(vals) <- features
  key <- data.frame(.u = long[[universe_col]], .v = long[[variant_col]],
                    vals, check.names = FALSE, stringsAsFactors = FALSE)

  # One row per (dataset, variant): the dataset's variant list. A variant that
  # appears in several regulatory elements keeps its largest feature value.
  safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  }
  uni <- dplyr::summarise(
    dplyr::group_by(key, dplyr::across(dplyr::all_of(c(".u", ".v")))),
    dplyr::across(dplyr::all_of(features), safe_max),
    .groups = "drop"
  )

  # Eq. 1: percentile within the dataset, p = #{u : x_u <= x_v} / |U|.
  pct <- function(x) {
    ok <- is.finite(x)
    out <- rep(NA_real_, length(x))
    if (any(ok)) out[ok] <- findInterval(x[ok], sort(x[ok])) / sum(ok)
    pmax(out, 1e-6)
  }
  uni <- dplyr::mutate(
    dplyr::group_by(uni, dplyr::across(dplyr::all_of(".u"))),
    dplyr::across(dplyr::all_of(features), pct, .names = "p_{.col}")
  )
  uni <- dplyr::ungroup(uni)

  # Eq. 2: geometric mean of the feature percentiles.
  pcols <- paste0("p_", features)
  uni$variant_score <- exp(rowMeans(log(as.matrix(uni[, pcols, drop = FALSE]))))

  size <- table(uni$.u)
  uni <- uni[, c(".u", ".v", pcols, "variant_score")]
  names(uni)[1:2] <- c(universe_col, variant_col)

  drop <- intersect(c(pcols, "variant_score"), names(long))
  if (length(drop)) long <- long[, setdiff(names(long), drop), drop = FALSE]
  out <- dplyr::left_join(dplyr::as_tibble(long), uni,
                          by = c(universe_col, variant_col))
  attr(out, "universe_size") <- size
  out
}

#' Aggregate variant scores to one raw score per sample
#'
#' Implements Eq. 3: the root mean square of the variant scores a sample
#' carries. The root mean square keeps the score on the same 0-1 scale as the
#' variant scores whatever the number of variants carried, so it does not
#' simply grow with variant burden, while still letting a few high-scoring
#' variants raise it.
#'
#' A variant seen more than once for the same sample (for example, a variant
#' that lies in two regulatory elements) is counted once.
#'
#' @param long Output of [pgvis_variant_score()].
#' @param by Grouping columns; `sample` plus whatever identifies the cohort
#'   and ancestry.
#' @param score_col Column holding the variant score.
#' @param variant_col Column identifying the variant.
#' @param metrics Optional metric columns whose per-sample mean is carried
#'   along (used by the plotting helpers); not part of the score.
#' @return One row per group with `n_variants` (variants contributing a score),
#'   `pgvis_raw` (the raw sample score S) and the requested metric means.
#' @export
pgvis_sample_scores <- function(long,
                                by = c("sample", "superpopulation", "cohort"),
                                score_col = "variant_score",
                                variant_col = "var_site",
                                metrics = c("delta_wi_int", "log_odds_ratio",
                                            "delta_af_alt")) {
  by <- intersect(by, names(long))
  if (!"sample" %in% by) stop("`by` must include 'sample'.", call. = FALSE)
  if (!score_col %in% names(long)) {
    stop("No column '", score_col, "'. Run pgvis_variant_score() first.",
         call. = FALSE)
  }
  metrics <- intersect(metrics, names(long))
  dedup <- unique(c(by, if (variant_col %in% names(long)) variant_col))
  d <- dplyr::distinct(dplyr::as_tibble(long),
                       dplyr::across(dplyr::all_of(dedup)), .keep_all = TRUE)
  out <- dplyr::summarise(
    dplyr::group_by(d, dplyr::across(dplyr::all_of(by))),
    n_variants = sum(!is.na(.data[[score_col]])),
    pgvis_raw  = sqrt(mean(.data[[score_col]]^2, na.rm = TRUE)),
    dplyr::across(dplyr::all_of(metrics), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )
  out$pgvis_raw[is.nan(out$pgvis_raw)] <- NA_real_
  out
}

#' Aggregate variant metrics to one row per sample (legacy)
#'
#' Averages raw metric columns per sample. Used by the legacy PCA composite
#' (`pgvis_score(method = "pca")`); the default score uses
#' [pgvis_variant_score()] and [pgvis_sample_scores()] instead.
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

#' Compute the PGViS score
#'
#' The default score (`method = "percentile"`, PGViS >= 0.2.0) standardises the
#' raw sample score `S` from [pgvis_sample_scores()] against the reference
#' cohort scored the same way (Eq. 4):
#'
#' \deqn{PGViS_i = (S_i - \mathrm{mean}_{ref}) / \mathrm{SD}_{ref}}{PGViS_i = (S_i - mean_ref) / SD_ref}
#'
#' The mean and standard deviation are taken over the reference samples of the
#' *same superpopulation* as sample `i`. Ancestry matching is not cosmetic:
#' allele frequency divergence from a reference panel is confounded with how
#' genetically distant a sample is from that panel, and control scores differ
#' by about one standard deviation across superpopulations.
#'
#' The full chain is
#'
#' ```
#' long   <- pgvis_variant_score(long)          # Eq. 1-2
#' sample <- pgvis_sample_scores(long)          # Eq. 3
#' scores <- pgvis_score(sample)                # Eq. 4
#' ```
#'
#' `method = "pca"` keeps the pre-0.2.0 composite: per-metric median/MAD
#' z-scores combined with PC1-derived weights. It is retained for
#' reproducibility only; its mean-aggregated metrics and prevalence weighting
#' follow the size of the variant list, and it includes `delta_af_alt`, which
#' is the variant-selection criterion.
#'
#' @param df One row per sample. For `method = "percentile"` it needs the
#'   `pgvis_raw` column from [pgvis_sample_scores()]; for `method = "pca"` the
#'   metric columns from [pgvis_aggregate_samples()].
#' @param metrics (`"pca"` only) Metric columns to combine.
#' @param reference_group Value of `group_col` identifying the reference
#'   cohort (e.g. `"1000Genomes"`).
#' @param group_col Column holding the cohort label.
#' @param strata Column defining the ancestry strata. Set to `NULL` to
#'   standardise against the whole reference cohort, which is only appropriate
#'   for an ancestry-homogeneous dataset.
#' @param weights (`"pca"` only) Optional named numeric vector to override the
#'   PCA-derived weights; see [pgvis_weights()].
#' @param method `"percentile"` (default) or `"pca"` (legacy).
#' @param reference (`"percentile"` only) Optional table of frozen reference
#'   statistics (columns `mu`, `sigma` and the `strata` column), e.g. the
#'   `"reference"` attribute of an earlier run. Use it to place a new cohort on
#'   the scale of a previous analysis without refitting.
#' @return One row per sample. For `"percentile"`: `pgvis_raw` and the
#'   standardised `pgvis`, with the reference statistics as attribute
#'   `"reference"`. For `"pca"`: the per-metric z-scores and `pgvis`, with the
#'   weights as attribute `"weights"`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   sample = paste0("S", 1:60),
#'   cohort = rep(c("1000Genomes", "tumour"), each = 30),
#'   superpopulation = "EUR", n_variants = 20L,
#'   pgvis_raw = c(rnorm(30, 0.5, 0.05), rnorm(30, 0.55, 0.05))
#' )
#' head(pgvis_score(d))
pgvis_score <- function(df,
                        metrics = c("delta_wi_int", "log_odds_ratio", "delta_af_alt"),
                        reference_group = "1000Genomes",
                        group_col = "cohort",
                        strata = "superpopulation",
                        weights = NULL,
                        method = c("percentile", "pca"),
                        reference = NULL) {
  method <- match.arg(method)
  if (method == "pca") {
    return(.pgvis_score_pca(df, metrics, reference_group, group_col, strata,
                            weights))
  }

  if (!"pgvis_raw" %in% names(df)) {
    stop("No 'pgvis_raw' column. Compute it with pgvis_variant_score() and ",
         "pgvis_sample_scores(), or use method = \"pca\" for the legacy ",
         "composite.", call. = FALSE)
  }
  if (!is.null(strata) && !strata %in% names(df)) {
    stop("No stratum column '", strata, "'. Pass strata = NULL to disable ancestry matching.",
         call. = FALSE)
  }

  if (is.null(reference)) {
    if (!group_col %in% names(df)) stop("No column '", group_col, "'.", call. = FALSE)
    if (!reference_group %in% df[[group_col]]) {
      stop("Reference group '", reference_group, "' not present in `", group_col,
           "`. The score is defined relative to a reference cohort; without one there is nothing to standardise against.",
           call. = FALSE)
    }
    ref <- dplyr::filter(df, .data[[group_col]] == reference_group)
    if (!is.null(strata)) ref <- dplyr::group_by(ref, .data[[strata]])
    ref <- dplyr::summarise(
      ref,
      mu    = mean(.data$pgvis_raw, na.rm = TRUE),
      sigma = stats::sd(.data$pgvis_raw, na.rm = TRUE),
      n_ref = sum(!is.na(.data$pgvis_raw)),
      .groups = "drop"
    )
  } else {
    need <- c("mu", "sigma", strata)
    miss <- setdiff(need, names(reference))
    if (length(miss)) stop("`reference` is missing: ", paste(miss, collapse = ", "),
                           call. = FALSE)
    ref <- reference
  }

  # A zero or undefined SD makes the z-score undefined; leave it NA rather than
  # invent a scale.
  bad <- !is.finite(ref$sigma) | ref$sigma == 0
  if (any(bad)) {
    .msg("Reference SD is undefined in ", sum(bad),
         " stratum/strata (too few reference samples or constant scores); ",
         "those samples score as NA.")
    ref$sigma[bad] <- NA_real_
  }

  ref_use <- ref[, intersect(c(strata, "mu", "sigma"), names(ref)), drop = FALSE]
  out <- if (is.null(strata)) {
    dplyr::bind_cols(dplyr::as_tibble(df), ref_use[rep(1L, nrow(df)), , drop = FALSE])
  } else {
    dplyr::left_join(dplyr::as_tibble(df), ref_use, by = strata)
  }
  out$pgvis <- (out$pgvis_raw - out$mu) / out$sigma

  keep <- c(intersect(c("sample", strata, group_col, "n_variants"), names(out)),
            "pgvis_raw", "pgvis")
  res <- out[, keep, drop = FALSE]
  attr(res, "reference") <- ref
  attr(res, "reference_group") <- reference_group
  attr(res, "method") <- "percentile"
  res
}

# Legacy composite (PGViS < 0.2.0): per-metric median/MAD z-scores combined with
# PC1-derived weights. Kept for reproducing earlier results.
.pgvis_score_pca <- function(df, metrics, reference_group, group_col, strata,
                             weights) {

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
  attr(res, "method") <- "pca"
  res
}

#' Derive metric weights from the first principal component
#'
#' Used by the legacy composite (`pgvis_score(method = "pca")`). Weights are
#' the absolute PC1 loadings, normalised to sum to one. Exported so that
#' weights fitted on a discovery cohort can be frozen and reused when scoring a
#' validation cohort, rather than being silently refitted.
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
