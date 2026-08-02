# ---------------------------------------------------------------------------
# Reference / alternate sequence construction.
#
# This is the step that turns a (regulatory region, variant) pair into the
# two sequences the model compares. It is strand-aware: the region sequence
# is extracted on the feature strand, the variant is applied in that frame,
# and both sequences are returned in the model's input orientation.
# ---------------------------------------------------------------------------

#' Apply a variant to a sequence
#'
#' Substitutes `alt` for `ref` at the expected offset. Because the region
#' start and the variant start come from separately produced BED files whose
#' coordinate conventions can disagree by a base, the reference allele is
#' searched for within a small window around the nominal offset and the
#' closest exact match is used.
#'
#' @param ref Reference allele.
#' @param alt Alternate allele. `""` or `"*"` performs a deletion.
#' @param seq The region sequence, in the same orientation as `ref`/`alt`.
#' @param pos_start,pos_end 1-based offsets of the variant within `seq`.
#' @param window Extra bases to search either side of `pos_start`. Defaults to
#'   the longer of the two alleles.
#' @return The mutated sequence, or `NA_character_` if `ref` could not be
#'   located (which is a data problem worth surfacing, not silently patching).
#' @export
#' @examples
#' mutate_sequence("A", "G", "TTTATTT", 4, 4)
#' mutate_sequence("AT", "", "TTTATTTT", 4, 5)
mutate_sequence <- function(ref, alt, seq, pos_start, pos_end, window = NULL) {
  if (is.na(seq) || is.na(ref)) return(NA_character_)
  alt <- if (is.na(alt)) "" else gsub("[*]", "", alt)
  ref_len <- nchar(ref)
  seq_len <- nchar(seq)
  if (ref_len == 0L || seq_len == 0L) return(NA_character_)

  window <- window %||% max(ref_len, nchar(alt), 1L)
  lo <- max(1L, pos_start - window)
  hi <- min(seq_len - ref_len + 1L, pos_end + window)
  if (hi < lo) return(NA_character_)

  offsets <- lo:hi
  matches <- substring(seq, offsets, offsets + ref_len - 1L) == ref
  if (!any(matches)) return(NA_character_)

  hit <- offsets[matches]
  best <- hit[which.min(abs(hit - pos_start))]
  paste0(substr(seq, 1L, best - 1L), alt,
         substr(seq, best + ref_len, seq_len))
}

#' Build reference and alternate sequences for scored variants
#'
#' @param variants A table with the region coordinates (`site_chr`,
#'   `site_start`, `site_end`, `site_strand`) and the variant coordinates
#'   (`var_chr`, `var_start`, `var_end`, `ref`, `alt`).
#' @param genome A `BSgenome` object. Defaults to
#'   `BSgenome.Hsapiens.UCSC.hg38`.
#' @param k K-mer size for the model input.
#' @param max_allele_len Drop variants with alleles longer than this.
#' @param drop_failures Remove rows where the alternate sequence could not be
#'   built. When `FALSE` they are returned with `NA` k-mers so you can inspect
#'   them.
#' @return `variants` with `ref_seq`, `alt_seq`, `ref_kmer` and `alt_kmer`
#'   columns added. Sequences are returned on the plus strand for TFBS and on
#'   the feature strand for splice sites, matching how the models were trained.
#' @export
pgvis_build_sequences <- function(variants,
                                  genome = NULL,
                                  k = 6L,
                                  max_allele_len = 10L,
                                  drop_failures = TRUE) {
  .need("GenomicRanges", "sequence extraction")
  .need("IRanges", "sequence extraction")
  .need("BSgenome", "sequence extraction")
  need_cols <- c("site_chr", "site_start", "site_end", "site_strand",
                 "var_start", "var_end", "ref", "alt")
  miss <- setdiff(need_cols, names(variants))
  if (length(miss)) {
    stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  if (is.null(genome)) {
    .need("BSgenome.Hsapiens.UCSC.hg38", "the default hg38 genome")
    genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38",
                               "BSgenome.Hsapiens.UCSC.hg38")
  }

  v <- dplyr::as_tibble(variants)
  v <- dplyr::filter(v,
                     nchar(.data$ref) <= max_allele_len,
                     nchar(gsub("[*]", "", .data$alt)) <= max_allele_len)
  # Strandless features are extracted on the plus strand.
  v$site_strand <- ifelse(is.na(v$site_strand) | v$site_strand %in% c(".", ""),
                          "+", v$site_strand)

  gr <- GenomicRanges::GRanges(
    seqnames = v$site_chr,
    ranges   = IRanges::IRanges(start = as.integer(v$site_start) + 1L,
                                end   = as.integer(v$site_end)),
    strand   = v$site_strand
  )
  v$ref_seq <- as.character(BSgenome::getSeq(genome, gr))

  minus <- v$site_strand == "-"
  L  <- nchar(v$ref_seq)
  vs <- as.integer(v$var_start)   # 0-based, BED
  ve <- as.integer(v$var_end)     # half-open, BED
  ss <- as.integer(v$site_start)  # 0-based, BED
  span <- pmax(1L, ve - vs)       # bases of REF covered

  # 1-based offset of the variant's first base within the strand-oriented
  # region sequence. On the minus strand the sequence is reverse
  # complemented, so the variant's last base becomes its first.
  v$mut_start <- ifelse(minus,
                        pmax(1L, L - (ve - ss) + 1L),
                        pmax(1L, vs - ss + 1L))
  v$mut_end <- pmin(L, v$mut_start + span - 1L)

  ref_o <- ifelse(minus, revcomp(v$ref), v$ref)
  alt_o <- ifelse(minus, revcomp(gsub("[*]", "", v$alt)), gsub("[*]", "", v$alt))

  v$alt_seq <- vapply(seq_len(nrow(v)), function(i) {
    mutate_sequence(ref_o[i], alt_o[i], v$ref_seq[i], v$mut_start[i], v$mut_end[i])
  }, character(1))

  v$ref_kmer <- seq2kmer(v$ref_seq, k)
  v$alt_kmer <- seq2kmer(v$alt_seq, k)

  n_fail <- sum(is.na(v$alt_kmer))
  if (n_fail) {
    .msg(sprintf("%d/%d variants: reference allele not found in the region sequence.",
                 n_fail, nrow(v)))
  }
  if (drop_failures) v <- dplyr::filter(v, !is.na(.data$ref_kmer), !is.na(.data$alt_kmer))
  v
}

#' Write DNABERT prediction input
#'
#' Writes one directory per feature containing `dev.tsv` (the two-column
#' k-mer/label file DNABERT's `run_finetune.py --do_predict` expects) and
#' `variants.tsv` (the full row metadata, needed to map predictions back).
#'
#' Rows alternate reference (label 0) then alternate (label 1), so that the
#' `pred_results.npy` array reshapes cleanly to two columns.
#'
#' @param variants Output of [pgvis_build_sequences()].
#' @param outdir Directory to write into, usually `cfg$dirs$seq`.
#' @param feature_col Column holding the feature/model name.
#' @return Character vector of the directories written, invisibly.
#' @export
pgvis_write_model_input <- function(variants, outdir, feature_col = "feature") {
  stopifnot(feature_col %in% names(variants))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  chunks <- split(variants, variants[[feature_col]])

  paths <- vapply(names(chunks), function(nm) {
    d <- file.path(outdir, nm)
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    df <- chunks[[nm]]
    utils::write.table(df, file.path(d, "variants.tsv"), sep = "\t",
                       row.names = FALSE, quote = FALSE)
    dev <- data.frame(
      Sequence = as.vector(rbind(df$ref_kmer, df$alt_kmer)),
      Label    = rep(c(0L, 1L), times = nrow(df))
    )
    utils::write.table(dev, file.path(d, "dev.tsv"), sep = "\t",
                       row.names = FALSE, quote = FALSE)
    d
  }, character(1), USE.NAMES = FALSE)

  .msg(sprintf("Wrote model input for %d feature(s) to %s", length(paths), outdir))
  invisible(paths)
}
