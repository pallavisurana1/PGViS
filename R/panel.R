# ---------------------------------------------------------------------------
# Reference panel support.
#
# The 1000 Genomes arm is the same pipeline pointed at a different VCF, but
# with three real differences, each handled here by an argument rather than a
# separate copy of the pipeline:
#
#   1. The panel VCF is already jointly genotyped, so stages 01-04 are skipped.
#   2. The intersection was run as `-a variants -b regions`, so the column
#      layout is reversed (`region_first = FALSE`).
#   3. Ancestry comes from the panel pedigree file rather than a clinical
#      table.
#
# All three are arguments here, so there is one pipeline, not two.
# ---------------------------------------------------------------------------

#' Configure a reference-panel project
#'
#' A [pgvis_config()] pre-set for a pre-genotyped panel such as 1000 Genomes:
#' the variant-calling stages are marked as not applicable and the cohort
#' label is fixed.
#'
#' @param root Project root.
#' @param vcf_dir Directory of per-chromosome panel VCFs.
#' @param label Cohort label used for these samples. This is the value you
#'   pass as `reference_group` to [pgvis_score()], so the two must agree.
#' @param ... Passed to [pgvis_config()].
#' @return A `pgvis_config` with `panel = TRUE`.
#' @export
#' @examples
#' cfg <- pgvis_panel_config(root = tempfile(), vcf_dir = ".")
#' cfg$cohorts
pgvis_panel_config <- function(root, vcf_dir, label = "1000Genomes", ...) {
  cfg <- pgvis_config(root = root, cohorts = label, ...)
  cfg$panel <- TRUE
  cfg$panel_vcf_dir <- vcf_dir
  cfg$dirs$joint <- vcf_dir   # stage 05 reads straight from the panel VCFs
  cfg
}

#' Load a 1000 Genomes-style pedigree / population file
#'
#' Reads the panel sample table and returns the two columns PGViS needs for
#' ancestry stratification. Column names are matched case-insensitively
#' because the published files are inconsistent about them
#' (`SampleID` / `sample` / `Sample`, `Superpopulation` / `super_pop`).
#'
#' @param path Path to the pedigree/population file.
#' @param sample_col,superpop_col Override the auto-detected column names.
#' @return A data frame with `sample`, `superpopulation` and, when present,
#'   `population` and `sex`.
#' @export
pgvis_load_pedigree <- function(path, sample_col = NULL, superpop_col = NULL) {
  .need("data.table", "reading the pedigree file")
  ped <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  nm <- tolower(names(ped))

  find_col <- function(given, candidates, what) {
    if (!is.null(given)) return(given)
    hit <- names(ped)[nm %in% candidates]
    if (!length(hit)) {
      stop("Could not find the ", what, " column. Looked for: ",
           paste(candidates, collapse = ", "),
           ". Found: ", paste(names(ped), collapse = ", "),
           ". Pass it explicitly.", call. = FALSE)
    }
    hit[[1]]
  }

  sc <- find_col(sample_col,
                 c("sampleid", "sample", "sample_id", "individualid", "iid"),
                 "sample ID")
  pc <- find_col(superpop_col,
                 c("superpopulation", "super_pop", "superpop", "super_population"),
                 "superpopulation")

  out <- data.frame(
    sample = as.character(ped[[sc]]),
    superpopulation = factor(toupper(trimws(as.character(ped[[pc]]))),
                             levels = c("AFR", "AMR", "EAS", "EUR", "SAS")),
    stringsAsFactors = FALSE
  )
  for (extra in c("population", "pop", "sex", "gender")) {
    hit <- names(ped)[nm == extra]
    if (length(hit)) out[[sub("gender", "sex", sub("^pop$", "population", extra))]] <-
        ped[[hit[[1]]]]
  }
  out <- unique(out)

  n_na <- sum(is.na(out$superpopulation))
  if (n_na) .msg(n_na, " sample(s) have an unrecognised superpopulation and will be excluded from scoring.")
  .msg(sprintf("Pedigree: %d samples across %d superpopulation(s).",
               nrow(out), length(stats::na.omit(unique(out$superpopulation)))))
  out
}

#' Harmonise a cohort table and a reference panel table for scoring
#'
#' Reduces both to the columns [pgvis_score()] needs, checks that the metric
#' columns and the ancestry levels actually overlap, and stacks them. The
#' failure mode this exists to prevent is scoring against a reference that
#' shares no superpopulation with the cohort, which produces all-`NA` z-scores
#' rather than an error.
#'
#' Each arm must already carry its raw sample score `pgvis_raw`, computed
#' against its *own* variant list with [pgvis_variant_score()] and
#' [pgvis_sample_scores()]. (Alternatively, stack the two arms' long variant
#' tables before scoring so a single call handles both.)
#'
#' @param cohort,panel Per-sample tables from [pgvis_sample_scores()] (or, for
#'   the legacy composite, [pgvis_aggregate_samples()]).
#' @param metrics Columns required in both. `NULL` (default) uses `pgvis_raw`
#'   when both tables have it, otherwise the three legacy metric columns.
#' @param panel_label Cohort label to assign to the panel rows.
#' @return A single stacked table ready for [pgvis_score()].
#' @export
pgvis_bind_panel <- function(cohort, panel,
                             metrics = NULL,
                             panel_label = "1000Genomes") {
  if (is.null(metrics)) {
    metrics <- if ("pgvis_raw" %in% names(cohort) && "pgvis_raw" %in% names(panel)) {
      "pgvis_raw"
    } else {
      c("delta_wi_int", "log_odds_ratio", "delta_af_alt")
    }
  }
  for (nm in c("cohort", "panel")) {
    d <- get(nm)
    miss <- setdiff(c("sample", metrics), names(d))
    if (length(miss)) {
      stop("`", nm, "` is missing: ", paste(miss, collapse = ", "), call. = FALSE)
    }
  }
  if (!"cohort" %in% names(panel)) panel$cohort <- panel_label

  if ("superpopulation" %in% names(cohort) && "superpopulation" %in% names(panel)) {
    shared <- intersect(
      as.character(stats::na.omit(unique(cohort$superpopulation))),
      as.character(stats::na.omit(unique(panel$superpopulation)))
    )
    if (!length(shared)) {
      stop("The cohort and the reference panel share no superpopulation. ",
           "Every z-score would be NA. Check the ancestry mapping on both sides.",
           call. = FALSE)
    }
    only_cohort <- setdiff(as.character(stats::na.omit(unique(cohort$superpopulation))), shared)
    if (length(only_cohort)) {
      .msg("No panel reference for superpopulation(s): ",
           paste(only_cohort, collapse = ", "),
           " -- those samples will score as NA.")
    }
  }

  keep <- unique(c("sample", "cohort", "superpopulation", "n_variants", metrics))
  bind <- function(d) {
    d <- as.data.frame(d)
    d[, intersect(keep, names(d)), drop = FALSE]
  }
  out <- dplyr::bind_rows(bind(cohort), bind(panel))
  .msg(sprintf("Stacked %d cohort + %d panel samples.",
               nrow(cohort), nrow(panel)))
  out
}
