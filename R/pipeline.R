# ---------------------------------------------------------------------------
# High-level workflow.
#
# Stages 1-6 are shell; stages 7-10 are R. The wrappers below give both a
# uniform interface so a whole project reads as a linear script.
# ---------------------------------------------------------------------------

#' Write the region manifest consumed by the intersection stage
#'
#' @param cfg A [pgvis_config()] whose `regions` slot is populated.
#' @param path Output path. Defaults to `<root>/regions.tsv`.
#' @param types Named character vector mapping feature name to
#'   `"splice"` or `"tfbs"`. Anything not named is inferred: `acceptor` and
#'   `donor` are splice, everything else is TFBS.
#' @return The manifest path, invisibly.
#' @export
pgvis_write_region_manifest <- function(cfg, path = NULL, types = NULL) {
  .assert_config(cfg)
  if (!length(cfg$regions)) {
    stop("cfg$regions is empty; nothing to intersect against.", call. = FALSE)
  }
  nm <- names(cfg$regions)
  ft <- ifelse(nm %in% c("acceptor", "donor"), "splice", "tfbs")
  if (!is.null(types)) {
    hit <- match(nm, names(types))
    ft[!is.na(hit)] <- types[hit[!is.na(hit)]]
  }
  path <- path %||% file.path(cfg$root, "regions.tsv")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data.frame(feature = nm, feature_type = ft,
               bed = unlist(cfg$regions, use.names = FALSE)),
    path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  .msg("Region manifest (", length(nm), " features) -> ", path)
  invisible(path)
}

#' Stage 1-6 wrappers
#'
#' Thin wrappers over [pgvis_render_script()] that fill in the variables each
#' shell stage needs from `cfg`. Each returns the path to a rendered,
#' executable script; pass `submit = TRUE` to run or queue it immediately.
#'
#' @param cfg A [pgvis_config()].
#' @param sample_sheet Headerless TSV: cohort, tissue, sample_id.
#' @param submit Run/queue the stage now.
#' @param array Optional SLURM array specification.
#' @param ... Extra shell variables, passed through to the template.
#' @name pgvis_shell_stages
#' @return Path to the rendered script, invisibly.
NULL

#' @rdname pgvis_shell_stages
#' @export
pgvis_merge_bams <- function(cfg, sample_sheet, submit = FALSE, array = NULL, ...) {
  pgvis_render_script("01_merge_bams.sh",
                      c(list(SAMPLE_SHEET = sample_sheet,
                             PGVIS_BASH_LIB = system.file("bash", package = "PGViS")),
                        list(...)),
                      cfg = cfg, submit = submit, array = array)
}

#' @rdname pgvis_shell_stages
#' @export
pgvis_call_gvcf <- function(cfg, sample_sheet, submit = FALSE,
                            array = NULL, ...) {
  if (is.null(cfg$reference_fasta) || is.na(cfg$reference_fasta)) {
    stop("cfg$reference_fasta must be set before calling variants.", call. = FALSE)
  }
  pgvis_render_script("02_call_gvcf.sh",
                      c(list(SAMPLE_SHEET = sample_sheet,
                             PGVIS_BASH_LIB = system.file("bash", package = "PGViS")),
                        list(...)),
                      cfg = cfg, submit = submit,
                      array = array %||% sprintf("1-%d", length(cfg$chromosomes)))
}

#' @rdname pgvis_shell_stages
#' @export
pgvis_combine_gvcf <- function(cfg, submit = FALSE, array = NULL, ...) {
  n <- length(cfg$cohorts) * length(cfg$chromosomes)
  pgvis_render_script("03_combine_gvcf.sh",
                      c(list(PGVIS_BASH_LIB = system.file("bash", package = "PGViS")),
                        list(...)),
                      cfg = cfg, submit = submit,
                      array = array %||% sprintf("1-%d", n))
}

#' @rdname pgvis_shell_stages
#' @export
pgvis_genotype_filter <- function(cfg, submit = FALSE, array = NULL, ...) {
  n <- length(cfg$cohorts) * length(cfg$chromosomes)
  pgvis_render_script("04_genotype_filter.sh",
                      c(list(PGVIS_BASH_LIB = system.file("bash", package = "PGViS")),
                        list(...)),
                      cfg = cfg, submit = submit,
                      array = array %||% sprintf("1-%d", n))
}

#' @rdname pgvis_shell_stages
#' @export
pgvis_vcf_to_bed <- function(cfg, submit = FALSE, array = NULL, ...) {
  pgvis_render_script("05_vcf2bed.sh",
                      c(list(PGVIS_BASH_LIB = system.file("bash", package = "PGViS")),
                        list(...)),
                      cfg = cfg, submit = submit, array = array)
}

#' @rdname pgvis_shell_stages
#' @param manifest Region manifest; created from `cfg` if `NULL`.
#' @export
pgvis_intersect <- function(cfg, manifest = NULL, submit = FALSE,
                            array = NULL, ...) {
  manifest <- manifest %||% pgvis_write_region_manifest(cfg)
  pgvis_render_script("06_intersect.sh",
                      c(list(REGION_MANIFEST = manifest,
                             PGVIS_BASH_LIB = system.file("bash", package = "PGViS")),
                        list(...)),
                      cfg = cfg, submit = submit,
                      array = array %||% sprintf("1-%d", length(cfg$regions)))
}

#' Load every intersected BED into one variant table
#'
#' Reads `<dbsnp dir>/<feature>/*_dbsnp.bed`, applies the schema, resolves
#' population frequencies and returns a single table.
#'
#' @param cfg A [pgvis_config()].
#' @param manifest Region manifest; read to recover each feature's type.
#' @param dir Directory of intersected BEDs. Defaults to `cfg$dirs$dbsnp`.
#' @param cohort_from_filename Regex with one capture group extracting the
#'   cohort label from the BED filename. Set to `NULL` to skip.
#' @param ... Passed to [pgvis_read_intersect()].
#' @return A `data.table` of annotated variants.
#' @export
pgvis_load_variants <- function(cfg, manifest = NULL, dir = NULL,
                                cohort_from_filename = "^([^_]+)_",
                                ...) {
  .assert_config(cfg)
  .need("data.table", "loading variants")
  dir <- dir %||% cfg$dirs$dbsnp
  manifest <- manifest %||% file.path(cfg$root, "regions.tsv")
  if (!file.exists(manifest)) {
    stop("Region manifest not found at ", manifest,
         "; run pgvis_write_region_manifest() first.", call. = FALSE)
  }
  man <- utils::read.delim(manifest, header = FALSE,
                           col.names = c("feature", "feature_type", "bed"),
                           stringsAsFactors = FALSE)

  out <- list()
  for (i in seq_len(nrow(man))) {
    fdir <- file.path(dir, man$feature[i])
    files <- list.files(fdir, pattern = "_dbsnp\\.bed$", full.names = TRUE)
    for (f in files) {
      df <- pgvis_read_intersect(f, feature_type = man$feature_type[i],
                                 feature = man$feature[i], ...)
      if (is.null(df) || !nrow(df)) next
      if (!is.null(cohort_from_filename)) {
        df$cohort <- stringr::str_match(basename(f), cohort_from_filename)[, 2]
      }
      df$source_file <- basename(f)
      out[[length(out) + 1L]] <- df
    }
  }
  if (!length(out)) stop("No non-empty intersected BEDs found under ", dir, call. = FALSE)

  res <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  res <- pgvis_resolve_population_af(res)
  .msg(sprintf("Loaded %d variant records across %d feature(s).",
               nrow(res), length(unique(res$feature))))
  data.table::as.data.table(res)
}

#' Run the R-side workflow end to end
#'
#' Sequence construction, model scoring, metrics, filtering, per-variant
#' percentile scoring, per-sample root-mean-square aggregation and
#' standardisation against the reference (see [pgvis_score()]), in one call. Every intermediate is
#' written to `cfg`'s directories, so you can stop after any stage and pick up
#' from the saved object.
#'
#' @param cfg A [pgvis_config()].
#' @param variants Output of [pgvis_load_variants()].
#' @param predictions Optional precomputed predictions. If `NULL`, PGViS
#'   builds the model input and attempts the Python backend; if that is
#'   unavailable it stops with instructions rather than guessing.
#' @param clinical Optional data frame with `sample` and either
#'   `superpopulation` or `race`, used for ancestry stratification.
#' @param reference_group Value of `cohort` naming the reference panel.
#' @param save Write intermediates to disk.
#' @return A list with `variants`, `metrics`, `per_sample` and `scores`.
#' @export
pgvis_run <- function(cfg, variants, predictions = NULL, clinical = NULL,
                      reference_group = "1000Genomes", save = TRUE) {
  .assert_config(cfg)
  pgvis_init(cfg)

  .msg("1/5 building reference and alternate sequences")
  seqs <- pgvis_build_sequences(variants,
                                max_allele_len = cfg$thresholds$max_indel_len)

  if (is.null(predictions)) {
    .msg("2/5 writing model input")
    dirs <- pgvis_write_model_input(seqs, cfg$dirs$seq)
    .msg("2/5 scoring with DeepVRegulome")
    predictions <- pgvis_predict_python(dirs, cfg$dirs$preds,
                                        python = cfg$python)
  } else {
    .msg("2/5 using supplied predictions")
  }

  .msg("3/5 computing per-variant metrics")
  metrics <- pgvis_variant_metrics(seqs, predictions, cfg = cfg)
  qc <- pgvis_qc(metrics, verbose = FALSE)
  if (length(qc$warnings)) {
    for (w in qc$warnings) warning("QC: ", w, call. = FALSE)
  }
  audit_before <- nrow(metrics)
  metrics <- pgvis_filter_disruptive(metrics, cfg = cfg)
  .msg("    ", nrow(metrics), " disruptive variant records retained of ", audit_before)
  if (!nrow(metrics)) {
    stop("No variants survived filtering. Loosen cfg$thresholds ",
         "(currently ref_prob > ", cfg$thresholds$ref_prob_min,
         ", alt_prob <= ", cfg$thresholds$alt_prob_max, ") ",
         "or check the QC report from pgvis_qc().", call. = FALSE)
  }
  if (save) saveRDS(metrics, file.path(cfg$dirs$metrics, "variant_metrics.rds"))

  .msg("4/5 expanding carriers and aggregating per sample")
  long <- pgvis_explode_carriers(metrics)
  if (!is.null(clinical)) {
    if (!"superpopulation" %in% names(clinical) && "race" %in% names(clinical)) {
      clinical$superpopulation <- pgvis_map_superpopulation(clinical$race)
    }
    long <- dplyr::left_join(long, clinical, by = "sample")
  }
  strata <- if ("superpopulation" %in% names(long)) "superpopulation" else NULL
  # Eq. 1-3: percentile within each cohort's own variant list, geometric mean
  # of |delta Wi| and LOR per variant, root mean square per sample.
  long <- pgvis_variant_score(long)
  per_sample <- pgvis_sample_scores(long, by = c("sample", strata, "cohort"))

  .msg("5/5 standardising against the ", reference_group, " reference")
  scores <- pgvis_score(per_sample, reference_group = reference_group,
                        strata = strata)   # Eq. 4
  if (save) {
    saveRDS(per_sample, file.path(cfg$dirs$score, "per_sample_metrics.rds"))
    saveRDS(scores, file.path(cfg$dirs$score, "pgvis_scores.rds"))
    utils::write.csv(scores, file.path(cfg$dirs$score, "pgvis_scores.csv"),
                     row.names = FALSE)
    pgvis_manifest(
      cfg, "pgvis_run",
      inputs  = list(variant_rows = nrow(variants),
                     sequence_rows = nrow(seqs),
                     prediction_rows = nrow(predictions),
                     retained_after_filter = nrow(metrics)),
      outputs = list(scores = file.path(cfg$dirs$score, "pgvis_scores.csv")),
      extra   = list(reference_group = reference_group,
                     strata = strata %||% "none",
                     score_method = "percentile / geometric mean / root mean square",
                     score_features = "abs(delta_wi_int), log_odds_ratio",
                     reference = paste(
                       sprintf("%s: mean=%.4f sd=%.4f",
                               if (is.null(strata)) "all" else
                                 as.character(attr(scores, "reference")[[strata]]),
                               attr(scores, "reference")$mu,
                               attr(scores, "reference")$sigma),
                       collapse = "; "),
                     n_scored_samples = nrow(scores),
                     qc_warnings = if (length(qc$warnings)) qc$warnings else "none")
    )
    .msg("Results written to ", cfg$dirs$score)
  }

  list(variants = seqs, metrics = metrics,
       per_sample = per_sample, scores = scores, qc = qc)
}
