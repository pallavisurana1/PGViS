# ---------------------------------------------------------------------------
# BED column schemas.
#
# The column names for every intersected BED are defined here once and
# assembled by `pgvis_schema()`. Because the awk emitter and these schemas are
# two halves of one contract, keeping the layout in a single place means a
# change to one has an obvious single counterpart to update.
# ---------------------------------------------------------------------------

.schema_parts <- list(

  # Emitted by inst/bash/vcf2bed_carriers.awk -- one row per ALT allele.
  cohort = c("var_chr", "var_start", "var_end", "ref", "alt", "allele_idx",
             "vtype", "ac", "af_alt", "an", "n_carrier", "carrier_samples"),

  # UCSC dbSnp155 BED (17 columns).
  dbsnp = c("dbsnp_chr", "dbsnp_start", "dbsnp_end", "rsid", "dbsnp_ref",
            "dbsnp_alt_count", "dbsnp_alt", "dbsnp_shift_bases",
            "dbsnp_freq_count", "dbsnp_freqs", "dbsnp_major", "dbsnp_minor",
            "dbsnp_func", "dbsnp_vtype", "dbsnp_flags", "dbsnp_raw_id",
            "dbsnp_score"),

  # GENCODE-derived splice donor/acceptor BED.
  splice = c("site_chr", "site_start", "site_end", "transcript", "site_type",
             "site_strand", "exon", "site_file"),

  # ENCODE-derived TFBS peak BED prepared for DNABERT (300 bp windows).
  tfbs = c("site_chr", "site_start", "site_end", "site_seq", "site_n",
           "site_file", "encode_id", "site_extra")
)

#' Column schema for an intersected BED file
#'
#' Returns the full column vector for a `bedtools intersect -wa -wb` product,
#' in the order the columns appear on disk.
#'
#' @param feature_type `"splice"` or `"tfbs"`.
#' @param with_dbsnp Has the file already been intersected against dbSNP?
#' @param region_first Was the regulatory-region BED passed as `-a`? The
#'   cohort pipeline uses `-a region -b variants`; the 1000 Genomes pipeline
#'   uses the reverse. Both layouts are supported so the two code paths can
#'   share one reader.
#' @return Character vector of column names.
#' @export
#' @examples
#' pgvis_schema("splice", with_dbsnp = TRUE)
pgvis_schema <- function(feature_type = c("splice", "tfbs"),
                         with_dbsnp = TRUE,
                         region_first = TRUE) {
  feature_type <- match.arg(feature_type)
  region <- .schema_parts[[feature_type]]
  variant <- .schema_parts$cohort
  parts <- if (region_first) list(region, variant) else list(variant, region)
  if (with_dbsnp) parts <- append(parts, list(.schema_parts$dbsnp),
                                  after = if (region_first) 2L else 1L)
  unlist(parts, use.names = FALSE)
}

#' Read an intersected BED into a tidy variant table
#'
#' Applies the schema, splits multi-allelic dbSNP records, harmonises allele
#' case, drops coordinate-discordant rows, and builds the canonical
#' `func_site` / `var_site` keys. This single function replaces the four
#' reader for every intersected BED, whatever feature type produced it.
#'
#' @param path Path to the BED file.
#' @param feature_type `"splice"` or `"tfbs"`.
#' @param feature Feature label (e.g. the TF name). Defaults to the name of
#'   the parent directory, which is how the intersection outputs are laid out.
#' @param with_dbsnp,region_first Passed to [pgvis_schema()].
#' @param require_allele_match If `TRUE`, keep only rows where the cohort REF
#'   and ALT agree with the dbSNP record. The cohort path historically matched
#'   on position only; the 1000 Genomes path matched on alleles too. Make the
#'   choice explicit rather than implicit in which script you ran.
#' @return A `data.table` of variants, or `NULL` for an empty/unreadable file.
#' @export
pgvis_read_intersect <- function(path,
                                 feature_type = c("splice", "tfbs"),
                                 feature = NULL,
                                 with_dbsnp = TRUE,
                                 region_first = TRUE,
                                 require_allele_match = FALSE) {
  feature_type <- match.arg(feature_type)
  .need("data.table", "reading intersected BED files")

  if (!file.exists(path) || file.size(path) == 0) return(NULL)
  df <- tryCatch(
    data.table::fread(path, header = FALSE, data.table = FALSE,
                      showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !nrow(df)) return(NULL)

  cols <- pgvis_schema(feature_type, with_dbsnp, region_first)
  if (ncol(df) != length(cols)) {
    stop(sprintf(
      "Column count mismatch in %s: file has %d columns, schema expects %d.\nDid the awk emitter or the region BED change? See ?pgvis_schema.",
      basename(path), ncol(df), length(cols)), call. = FALSE)
  }
  names(df) <- cols

  feature <- feature %||% basename(dirname(path))
  df$feature <- feature
  df$feature_type <- feature_type
  if (feature_type == "tfbs") df$site_strand <- "."

  df <- dplyr::mutate(
    df,
    ref        = toupper(trimws(.data$ref)),
    alt        = toupper(trimws(.data$alt)),
    af_alt     = suppressWarnings(as.numeric(ifelse(.data$af_alt == ".", NA, .data$af_alt))),
    af_ref     = 1 - .data$af_alt,
    n_carrier  = suppressWarnings(as.integer(.data$n_carrier))
  )

  if (with_dbsnp) {
    df <- df |>
      dplyr::mutate(
        dbsnp_alt = toupper(trimws(sub(",$", "", .data$dbsnp_alt))),
        dbsnp_ref = toupper(trimws(.data$dbsnp_ref))
      ) |>
      tidyr::separate_rows("dbsnp_alt", sep = ",") |>
      dplyr::filter(.data$var_start == .data$dbsnp_start)
    if (require_allele_match) {
      df <- dplyr::filter(df, .data$ref == .data$dbsnp_ref,
                          .data$alt == .data$dbsnp_alt)
    }
  }

  df$func_site <- func_site_key(df$site_chr, df$site_start, df$site_end,
                                df$site_strand,
                                stranded = feature_type == "splice")
  df$var_site  <- variant_key(df$var_chr, df$var_start, df$var_end,
                              df$ref, df$alt)

  data.table::as.data.table(df)
}
