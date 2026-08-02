# ---------------------------------------------------------------------------
# Population allele frequency resolution.
#
# UCSC's dbSnp155 track packs per-project frequencies into comma-delimited
# strings whose element order is fixed by the track schema. The position of
# each project (TOPMED = 3, GnoMAD = 13, ...) and the order in which they are
# consulted are both arguments here, so the resolution logic is defined once
# and stays configurable.
# ---------------------------------------------------------------------------

#' Default dbSNP frequency-source index table
#'
#' Positions of each project inside the comma-delimited
#' `minorAlleleFreq` / `majorAllele` / `minorAllele` fields of the UCSC
#' dbSnp155 BED track.
#'
#' @return A named integer vector.
#' @export
#' @examples
#' pgvis_freq_sources()
pgvis_freq_sources <- function() {
  c(TOPMED = 3L, GnoMAD = 13L, `1000G` = 1L, dbGAP = 2L)
}

.nth_csv <- function(x, i) {
  vapply(strsplit(as.character(x), ",", fixed = TRUE),
         function(v) if (length(v) >= i) v[[i]] else NA_character_,
         character(1), USE.NAMES = FALSE)
}

#' Resolve reference and alternate population frequencies
#'
#' Walks the frequency sources in priority order and takes the first one whose
#' reported minor allele matches either the reference or the alternate allele
#' of the variant. The matched minor-allele frequency is assigned to whichever
#' allele it corresponds to, and `1 - maf` to the other.
#'
#' @param df A table produced by [pgvis_read_intersect()], containing
#'   `dbsnp_freqs`, `dbsnp_major`, `dbsnp_minor`, `dbsnp_ref` and `dbsnp_alt`.
#' @param sources Named integer vector of column indices; see
#'   [pgvis_freq_sources()].
#' @param priority Order in which sources are consulted. Defaults to the order
#'   of `sources`.
#' @param keep_raw Keep the per-source `maf_*` / `minor_*` columns for
#'   auditing.
#' @return `df` with three new columns: `pop_af_alt`, `pop_af_ref` and
#'   `freq_source` (the name of the source used, or `"none"`).
#' @export
#' @examples
#' toy <- data.frame(
#'   dbsnp_freqs = "0.1,0.2,0.31,,,,,,,,,,0.30",
#'   dbsnp_minor = "A,A,A,,,,,,,,,,A",
#'   dbsnp_major = "G,G,G,,,,,,,,,,G",
#'   dbsnp_ref = "G", dbsnp_alt = "A"
#' )
#' pgvis_resolve_population_af(toy)[, c("pop_af_alt", "freq_source")]
pgvis_resolve_population_af <- function(df,
                                        sources = pgvis_freq_sources(),
                                        priority = names(sources),
                                        keep_raw = FALSE) {
  stopifnot(all(c("dbsnp_freqs", "dbsnp_minor", "dbsnp_ref", "dbsnp_alt")
                %in% names(df)))
  priority <- intersect(priority, names(sources))
  if (!length(priority)) stop("No usable frequency sources.", call. = FALSE)

  n <- nrow(df)
  af_alt <- rep(NA_real_, n)
  af_ref <- rep(NA_real_, n)
  src    <- rep("none", n)

  ref <- toupper(trimws(df$dbsnp_ref))
  alt <- toupper(trimws(df$dbsnp_alt))

  for (s in priority) {
    i <- sources[[s]]
    maf   <- suppressWarnings(as.numeric(.nth_csv(df$dbsnp_freqs, i)))
    maf[!is.finite(maf)] <- NA_real_
    minor <- toupper(trimws(.nth_csv(df$dbsnp_minor, i)))

    if (keep_raw) {
      df[[paste0("maf_", s)]]   <- maf
      df[[paste0("minor_", s)]] <- minor
    }

    unfilled <- src == "none" & !is.na(maf)
    hit_alt <- unfilled & !is.na(minor) & minor == alt
    hit_ref <- unfilled & !is.na(minor) & minor == ref

    af_alt[hit_alt] <- maf[hit_alt]
    af_ref[hit_alt] <- 1 - maf[hit_alt]
    af_alt[hit_ref] <- 1 - maf[hit_ref]
    af_ref[hit_ref] <- maf[hit_ref]
    src[hit_alt | hit_ref] <- s
  }

  df$pop_af_alt  <- af_alt
  df$pop_af_ref  <- af_ref
  df$freq_source <- src
  df
}
