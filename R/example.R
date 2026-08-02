# ---------------------------------------------------------------------------
# Bundled example data.
#
# A small synthetic fixture so the vignette, the examples and CI can exercise
# the whole R-side workflow without a reference genome, a dbSNP track, or a
# GPU. The numbers are not biologically meaningful; they exist so the code
# paths are covered and so you can see the expected shape of each table.
# ---------------------------------------------------------------------------

#' Load a bundled example dataset
#'
#' @param which One of:
#'   \describe{
#'     \item{`"variants"`}{An annotated variant table of the kind
#'       [pgvis_load_variants()] returns: 240 rows over 4 features, with
#'       cohort and population allele frequencies and carrier lists.}
#'     \item{`"predictions"`}{Matching model probabilities, joinable on
#'       `feature`/`func_site`/`var_site`.}
#'     \item{`"clinical"`}{Per-sample race and survival follow-up.}
#'     \item{`"pedigree"`}{A 1000 Genomes-style population file.}
#'     \item{`"samples"`}{A three-column sample sheet.}
#'     \item{`"regions"`}{A region manifest.}
#'   }
#' @return A data frame.
#' @export
#' @examples
#' v <- pgvis_example("variants")
#' str(v[, 1:8])
#' table(v$feature)
pgvis_example <- function(which = c("variants", "predictions", "clinical",
                                    "pedigree", "samples", "regions")) {
  which <- match.arg(which)
  f <- switch(which,
              variants    = "toy_variants.tsv",
              predictions = "toy_predictions.tsv",
              clinical    = "toy_clinical.csv",
              pedigree    = "toy_pedigree.txt",
              samples     = "toy_samples.tsv",
              regions     = "regions.tsv")
  path <- system.file("extdata", f, package = "PGViS")
  if (!nzchar(path)) stop("Example file '", f, "' not found. Reinstall PGViS.",
                          call. = FALSE)

  if (which == "clinical") {
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  }
  if (which == "samples") {
    return(utils::read.delim(path, header = FALSE,
                             col.names = c("cohort", "tissue", "sample"),
                             stringsAsFactors = FALSE))
  }
  if (which == "regions") {
    return(utils::read.delim(path, header = FALSE,
                             col.names = c("feature", "feature_type", "bed"),
                             stringsAsFactors = FALSE))
  }
  utils::read.delim(path, stringsAsFactors = FALSE)
}

#' Run the whole R-side workflow on the bundled example
#'
#' A single call that goes from the example variant table to scored samples,
#' skipping the genome and the model. Useful as a smoke test after installing,
#' and as a template to read alongside your own analysis.
#'
#' @param verbose Print progress.
#' @return A list with `metrics`, `per_sample` and `scores`.
#' @export
#' @examples
#' res <- pgvis_example_run()
#' head(res$scores)
#' attr(res$scores, "weights")
pgvis_example_run <- function(verbose = TRUE) {
  old <- getOption("PGViS.quiet")
  options(PGViS.quiet = !verbose)
  on.exit(options(PGViS.quiet = old), add = TRUE)

  variants <- pgvis_example("variants")
  preds    <- pgvis_example("predictions")
  clin     <- pgvis_example("clinical")
  ped      <- pgvis_load_pedigree(
    system.file("extdata", "toy_pedigree.txt", package = "PGViS"))

  metrics <- pgvis_variant_metrics(variants, preds)
  metrics <- pgvis_filter_disruptive(metrics,
                                     ref_prob_min = 0.5, alt_prob_max = 0.5,
                                     min_carrier_frac = 0)

  long <- pgvis_explode_carriers(metrics)
  clin$superpopulation <- pgvis_map_superpopulation(clin$race)
  long <- dplyr::left_join(long, clin[, c("sample", "superpopulation")],
                           by = "sample")

  cohort_ps <- pgvis_aggregate_samples(long)
  cohort_ps <- cohort_ps[!is.na(cohort_ps$superpopulation), , drop = FALSE]

  # Synthesise a reference arm from the pedigree so the example can be scored
  # end to end. In a real project this comes from running the same pipeline
  # over the panel VCFs; see pgvis_panel_config().
  set.seed(11)
  panel_ps <- data.frame(
    sample          = ped$sample,
    superpopulation = ped$superpopulation,
    cohort          = "1000Genomes",
    n_variants      = 10L,
    delta_wi_int    = stats::rnorm(nrow(ped)),
    log_odds_ratio  = stats::rnorm(nrow(ped)),
    delta_af_alt    = stats::rnorm(nrow(ped)),
    stringsAsFactors = FALSE
  )

  stacked <- pgvis_bind_panel(cohort_ps, panel_ps)
  scores  <- pgvis_score(stacked, reference_group = "1000Genomes")

  list(metrics = metrics, per_sample = stacked, scores = scores)
}
