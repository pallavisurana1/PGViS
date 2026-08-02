# ---------------------------------------------------------------------------
# Cohort bookkeeping.
#
# Cohort identity often travels as a composite string like "tcga_luad_smoke".
# Deriving study, cancer type and stratum from it ad hoc invites subtle bugs:
# two slightly different regexes will split a group differently and nothing
# complains. Parsing happens once, here, and the components are kept as real
# columns from then on.
# ---------------------------------------------------------------------------

#' Build or validate a sample sheet
#'
#' The shell stages consume a headerless three-column TSV: cohort, tissue,
#' sample_id. This writes one and checks it is usable before you queue a
#' thousand HaplotypeCaller jobs behind it.
#'
#' @param x A data frame with columns `cohort`, `tissue`, `sample`, or a path
#'   to an existing sheet to validate.
#' @param path Where to write. `NULL` validates without writing.
#' @param bam_dir If supplied, check that a BAM exists for every row and warn
#'   about the ones that do not.
#' @return The validated data frame, invisibly.
#' @export
#' @examples
#' pgvis_sample_sheet(data.frame(cohort = "luad", tissue = "blood",
#'                               sample = "S1"))
pgvis_sample_sheet <- function(x, path = NULL, bam_dir = NULL) {
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    x <- utils::read.delim(x, header = FALSE,
                           col.names = c("cohort", "tissue", "sample"),
                           stringsAsFactors = FALSE)
  }
  need <- c("cohort", "tissue", "sample")
  miss <- setdiff(need, names(x))
  if (length(miss)) {
    stop("Sample sheet needs columns: ", paste(need, collapse = ", "),
         ". Missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  x <- x[, need, drop = FALSE]

  # A tab inside a field silently shifts every downstream column.
  bad <- vapply(x, function(col) any(grepl("[\t\r\n]", col)), logical(1))
  if (any(bad)) {
    stop("Whitespace that would break the TSV in column(s): ",
         paste(names(bad)[bad], collapse = ", "), call. = FALSE)
  }
  dup <- duplicated(x[, c("cohort", "tissue", "sample")])
  if (any(dup)) {
    stop(sum(dup), " duplicated cohort/tissue/sample row(s); ",
         "the GVCF path is built from those three fields, so duplicates would overwrite each other.",
         call. = FALSE)
  }

  if (!is.null(bam_dir)) {
    found <- vapply(seq_len(nrow(x)), function(i) {
      d <- file.path(bam_dir, x$cohort[i], x$tissue[i], x$sample[i])
      length(list.files(d, pattern = "\\.bam$", recursive = TRUE)) > 0
    }, logical(1))
    if (any(!found)) {
      warning(sum(!found), "/", nrow(x), " samples have no BAM under ",
              bam_dir, ". First few: ",
              paste(utils::head(x$sample[!found], 3), collapse = ", "),
              call. = FALSE)
    }
  }

  if (!is.null(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(x, path, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
    .msg("Sample sheet (", nrow(x), " rows) -> ", path)
  }
  invisible(x)
}

#' Split a composite source label into its components
#'
#' Turns `"tcga_luad_smoke"` into study `tcga`, cancer `luad`, stratum
#' `smoke`. Unrecognised components become `NA` rather than being pattern-
#' matched into the wrong slot.
#'
#' @param label Character vector of composite labels.
#' @param studies,cancers,strata Recognised values for each position.
#' @return A data frame with `study`, `cancer`, `stratum`.
#' @export
#' @examples
#' pgvis_parse_label(c("tcga_luad_smoke", "cptac_lusc", "1000Genomes"))
pgvis_parse_label <- function(label,
                              studies = c("tcga", "cptac"),
                              cancers = c("luad", "lusc"),
                              strata  = c("smoke", "nonsmoke")) {
  lab <- tolower(as.character(label))
  pick <- function(vals) {
    vapply(lab, function(s) {
      hit <- vals[vapply(vals, function(v) grepl(v, s, fixed = TRUE), logical(1))]
      if (length(hit)) hit[[1]] else NA_character_
    }, character(1), USE.NAMES = FALSE)
  }
  data.frame(study = pick(studies), cancer = pick(cancers),
             stratum = pick(strata), stringsAsFactors = FALSE)
}

#' Combine per-cohort result files into one table
#'
#' Reads a set of CSV/TSV/parquet files and stamps each row with where it came
#' from. Replaces the `map_dfr()` + inline `str_extract()` block that was
#' rewritten in `5B_LOD_filt.R`, `3_5A_analysis_Smoke.R` and elsewhere.
#'
#' @param files Character vector of paths.
#' @param labels Source label per file. Defaults to the file's parent
#'   directory path components joined by `_`, which is usually what you want
#'   given the `study/cancer/stratum/` layout.
#' @param parse_label Also split the label into `study`/`cancer`/`stratum`.
#' @param ... Passed to [pgvis_parse_label()].
#' @return A single `data.table` with `source_label` and `source_file`.
#' @export
pgvis_combine_results <- function(files, labels = NULL, parse_label = TRUE, ...) {
  .need("data.table", "combining result files")
  files <- files[file.exists(files)]
  if (!length(files)) stop("None of the supplied files exist.", call. = FALSE)

  labels <- labels %||% vapply(files, function(f) {
    p <- strsplit(normalizePath(f, mustWork = FALSE), .Platform$file.sep)[[1]]
    paste(utils::tail(p[-length(p)], 3), collapse = "_")
  }, character(1), USE.NAMES = FALSE)

  out <- lapply(seq_along(files), function(i) {
    f <- files[i]
    df <- if (grepl("\\.parquet$", f)) {
      .need("arrow", "reading parquet")
      as.data.frame(arrow::read_parquet(f))
    } else {
      data.table::fread(f, data.table = FALSE, showProgress = FALSE)
    }
    if (!nrow(df)) return(NULL)
    df$source_label <- labels[i]
    df$source_file  <- f
    df
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) stop("Every file was empty.", call. = FALSE)

  res <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  if (parse_label) res <- cbind(res, pgvis_parse_label(res$source_label, ...))
  .msg(sprintf("Combined %d file(s), %d rows, %d distinct source label(s).",
               length(out), nrow(res), length(unique(res$source_label))))
  res[]
}

#' Keep only variants replicated across independent studies
#'
#' A variant seen in TCGA and again in CPTAC is far more credible than one
#' seen in either alone. This keeps sites observed in at least `min_studies`
#' distinct studies, within each stratum.
#'
#' @param df A combined result table with a `study` column.
#' @param keys Columns identifying a variant.
#' @param within Additional grouping columns; replication is required
#'   separately within each level (so a smoking-stratified analysis does not
#'   count a smoker hit and a non-smoker hit as replication).
#' @param study_col Column holding the study identifier.
#' @param min_studies Minimum distinct studies.
#' @return The filtered table, with a `n_studies` column.
#' @export
pgvis_filter_replicated <- function(df,
                                    keys = c("func_site", "var_site"),
                                    within = "cancer",
                                    study_col = "study",
                                    min_studies = 2L) {
  keys   <- intersect(keys, names(df))
  within <- intersect(within, names(df))
  if (!study_col %in% names(df)) {
    stop("No `", study_col, "` column. Run pgvis_combine_results() first, ",
         "or pass study_col.", call. = FALSE)
  }
  if (!length(keys)) stop("None of the key columns are present.", call. = FALSE)

  before <- nrow(df)
  out <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(within, keys)))) |>
    dplyr::mutate(n_studies = dplyr::n_distinct(.data[[study_col]])) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$n_studies >= min_studies)

  .msg(sprintf("Replication filter (>= %d studies): %d -> %d rows (%d unique sites kept).",
               min_studies, before, nrow(out),
               nrow(unique(out[, keys, drop = FALSE]))))
  out
}

#' Summarise a variant table by group
#'
#' Counts unique functional sites, variant sites and rsIDs per group. The
#' table you want in a methods section.
#'
#' @param df A variant table.
#' @param by Grouping columns.
#' @return A `data.frame`, one row per group.
#' @export
pgvis_summarise_by <- function(df, by = c("study", "cancer", "feature_type")) {
  by <- intersect(by, names(df))
  if (!length(by)) return(pgvis_count_sites(df, "all"))
  chunks <- split(df, interaction(df[, by, drop = FALSE], drop = TRUE))
  res <- do.call(rbind, lapply(names(chunks), function(nm) {
    r <- pgvis_count_sites(chunks[[nm]], nm)
    cbind(chunks[[nm]][1, by, drop = FALSE], r[, -1, drop = FALSE])
  }))
  rownames(res) <- NULL
  res
}
