# ---------------------------------------------------------------------------
# Provenance and QC.
#
# "Which thresholds produced this table?" is a question that gets asked months
# after the run, when nobody remembers. A run manifest is written alongside
# every result so the settings travel with the data rather than living only in
# whatever code happened to be current at the time.
# ---------------------------------------------------------------------------

#' Record the settings and inputs behind a result
#'
#' @param cfg A [pgvis_config()].
#' @param stage Short stage name.
#' @param inputs Named list of input descriptions (paths, row counts).
#' @param outputs Named list of output paths.
#' @param extra Anything else worth recording.
#' @param path Where to write the JSON-ish manifest. Defaults to the logs
#'   directory, timestamped.
#' @return The manifest as a list, invisibly.
#' @export
pgvis_manifest <- function(cfg, stage, inputs = list(), outputs = list(),
                           extra = list(), path = NULL) {
  .assert_config(cfg)
  man <- list(
    stage      = stage,
    timestamp  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    pgvis      = as.character(utils::packageVersion("PGViS")),
    r_version  = paste(R.version$major, R.version$minor, sep = "."),
    genome     = cfg$genome,
    cohorts    = cfg$cohorts,
    thresholds = cfg$thresholds,
    n_features = length(cfg$regions),
    inputs     = inputs,
    outputs    = outputs,
    extra      = extra
  )
  path <- path %||% file.path(
    cfg$dirs$logs,
    sprintf("manifest_%s_%s.txt", stage, format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(.flatten_manifest(man), path)
  .msg("Manifest -> ", path)
  invisible(man)
}

.flatten_manifest <- function(x, prefix = "") {
  out <- character(0)
  for (nm in names(x)) {
    v <- x[[nm]]
    key <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
    if (is.list(v)) {
      out <- c(out, .flatten_manifest(v, key))
    } else {
      out <- c(out, sprintf("%s: %s", key,
                            paste(format(v, trim = TRUE), collapse = ", ")))
    }
  }
  out
}

#' Quality-control summary of a variant table
#'
#' Reports the things that most often turn out to be wrong: features with no
#' variants at all, missing population frequencies, model probabilities that
#' never move, and carrier counts that look like the genotype columns were
#' misparsed.
#'
#' @param df A variant or metric table.
#' @param verbose Print the report as well as returning it.
#' @return A named list of QC findings, invisibly when `verbose`.
#' @export
pgvis_qc <- function(df, verbose = TRUE) {
  n <- nrow(df)
  qc <- list(n_rows = n)

  has <- function(x) x %in% names(df)
  pct <- function(k) sprintf("%.1f%%", 100 * k / max(n, 1))

  if (has("feature")) {
    tab <- table(df$feature)
    qc$n_features <- length(tab)
    qc$sparsest_features <- utils::head(sort(tab), 5)
  }
  if (has("freq_source")) {
    qc$freq_source <- table(df$freq_source)
    none <- sum(df$freq_source == "none")
    qc$no_population_af <- none
    qc$no_population_af_pct <- pct(none)
  }
  if (has("pop_af_alt")) qc$missing_pop_af <- sum(is.na(df$pop_af_alt))
  if (has("af_alt"))     qc$missing_cohort_af <- sum(is.na(df$af_alt))

  if (has("p_ref") && has("p_alt")) {
    delta <- abs(df$p_ref - df$p_alt)
    qc$median_prob_shift <- stats::median(delta, na.rm = TRUE)
    qc$unmoved_predictions <- sum(delta < 1e-9, na.rm = TRUE)
    # Every probability identical usually means the .npy was reshaped wrongly.
    qc$degenerate_predictions <-
      length(unique(round(c(df$p_ref, df$p_alt), 6))) <= 2L
  }
  if (has("n_carrier")) {
    qc$carrier_range <- range(df$n_carrier, na.rm = TRUE)
    qc$zero_carrier_variants <- sum(df$n_carrier == 0, na.rm = TRUE)
  }
  if (has("delta_wi_int")) {
    qc$delta_wi_na <- sum(is.na(df$delta_wi_int))
  }

  warnings <- character(0)
  if (!is.null(qc$no_population_af) && qc$no_population_af > 0.5 * n) {
    warnings <- c(warnings,
      "Over half the variants have no matched population frequency. Check that the dbSNP build matches the genome and that the frequency-source indices are right for your track version (see pgvis_freq_sources()).")
  }
  if (isTRUE(qc$degenerate_predictions)) {
    warnings <- c(warnings,
      "Model probabilities take almost no distinct values. The prediction array was probably reshaped incorrectly, or the wrong model was applied.")
  }
  if (!is.null(qc$zero_carrier_variants) && qc$zero_carrier_variants > 0) {
    warnings <- c(warnings, sprintf(
      "%d variant(s) have zero carriers -- they cannot contribute to any sample score and are probably a genotype-parsing artefact.",
      qc$zero_carrier_variants))
  }
  if (!is.null(qc$unmoved_predictions) && qc$unmoved_predictions > 0.9 * n) {
    warnings <- c(warnings,
      "Almost no variant changes the model's probability. Check that alt_seq really differs from ref_seq (pgvis_build_sequences reports failures).")
  }
  qc$warnings <- warnings

  if (verbose) {
    cat("PGViS QC\n--------\n")
    cat("rows:", qc$n_rows, "\n")
    for (nm in setdiff(names(qc), c("warnings", "n_rows"))) {
      v <- qc[[nm]]
      cat(sprintf("%-24s %s\n", paste0(nm, ":"),
                  paste(utils::capture.output(print(v)), collapse = " ")))
    }
    if (length(warnings)) {
      cat("\nWarnings\n--------\n")
      for (w in warnings) cat("* ", strwrap(w, exdent = 2), "\n", sep = "")
    } else {
      cat("\nNo QC warnings.\n")
    }
    return(invisible(qc))
  }
  qc
}
