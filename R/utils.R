#' @keywords internal
"_PACKAGE"

`%||%` <- function(x, y) if (is.null(x)) y else x

.msg <- function(...) {
  if (isTRUE(getOption("PGViS.quiet", FALSE))) return(invisible(NULL))
  message("[PGViS] ", ...)
}

.need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s. Install it with install.packages('%s').",
                 pkg, what, pkg), call. = FALSE)
  }
  invisible(TRUE)
}

.need_bin <- function(bin) {
  if (nzchar(Sys.which(bin)[[1]])) return(invisible(TRUE))
  stop(sprintf("External tool '%s' was not found on PATH.", bin), call. = FALSE)
}

#' Symmetric difference of two vectors
#'
#' Convenience helper: elements in exactly one of two sets.
#'
#' @param x,y Vectors.
#' @return Elements present in exactly one of `x`, `y`.
#' @export
outersect <- function(x, y) setdiff(union(x, y), intersect(x, y))

# ---------------------------------------------------------------------------
# Coordinate keys
#
# Coordinate keys are the join keys between the variant table and the model
# predictions, so they have to match exactly on both sides. Building them
# inline invites small spelling differences that break joins silently. There
# is exactly one constructor for each kind of key, used everywhere.
# ---------------------------------------------------------------------------

#' Build a functional-site key
#'
#' Produces `chr:start-end(strand)`. Strandless features (TFBS) get an empty
#' `()` so that keys are comparable across feature types.
#'
#' @param chr,start,end Coordinate vectors.
#' @param strand Strand vector; `"."`, `NA` and `""` are all treated as
#'   strandless.
#' @param stranded Logical vector (recycled): should the strand be kept?
#'   Splice donor/acceptor sites are stranded, TFBS are not.
#' @return Character vector of site keys.
#' @export
#' @examples
#' func_site_key("chr1", 100, 180, "+")
#' func_site_key("chr1", 100, 180, "+", stranded = FALSE)
func_site_key <- function(chr, start, end, strand = ".", stranded = TRUE) {
  strand <- as.character(strand)
  strand[is.na(strand) | strand %in% c(".", "")] <- ""
  strand[!stranded] <- ""
  paste0(chr, ":", start, "-", end, "(", strand, ")")
}

#' Build a variant key
#'
#' Produces `chr:start-end(REF>ALT)`. Empty, `NA` or non-ACGT alternate
#' alleles are normalised to `*`, matching the convention used by the DNABERT
#' prediction tables.
#'
#' @param chr,start,end Coordinate vectors.
#' @param ref,alt Allele vectors.
#' @return Character vector of variant keys.
#' @export
#' @examples
#' variant_key("chr1", 120, 121, "A", "G")
#' variant_key("chr1", 120, 121, "A", NA)
variant_key <- function(chr, start, end, ref, alt) {
  alt <- toupper(as.character(alt))
  alt[is.na(alt) | !grepl("^[ACGTN]+$", alt)] <- "*"
  paste0(chr, ":", start, "-", end, "(", toupper(ref), ">", alt, ")")
}

#' Split a site or variant key back into columns
#'
#' @param key Character vector produced by [func_site_key()] or
#'   [variant_key()].
#' @param type `"site"` or `"variant"`.
#' @return A `data.frame`.
#' @export
parse_key <- function(key, type = c("site", "variant")) {
  type <- match.arg(type)
  chr   <- stringr::str_extract(key, "^[^:]+")
  start <- as.integer(stringr::str_extract(key, "(?<=:)\\d+"))
  end   <- as.integer(stringr::str_extract(key, "(?<=-)\\d+"))
  inner <- stringr::str_extract(key, "(?<=\\()[^)]*(?=\\))")
  if (type == "site") {
    strand <- ifelse(is.na(inner) | inner == "", ".", inner)
    data.frame(chr = chr, start = start, end = end, strand = strand,
               stringsAsFactors = FALSE)
  } else {
    data.frame(chr = chr, start = start, end = end,
               ref = sub(">.*$", "", inner),
               alt = sub("^.*>", "", inner),
               stringsAsFactors = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

#' Reverse complement a DNA string
#'
#' A dependency-free implementation, so that the metric and sequence layers of
#' PGViS do not force a Bioconductor install. Handles ambiguity codes.
#'
#' @param seq Character vector of DNA sequences.
#' @return Character vector of the same length.
#' @export
#' @examples
#' revcomp("ACCGGT")
revcomp <- function(seq) {
  map <- c(A = "T", C = "G", G = "C", T = "A", N = "N",
           R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
           B = "V", V = "B", D = "H", H = "D")
  vapply(seq, function(s) {
    if (is.na(s)) return(NA_character_)
    ch <- rev(strsplit(toupper(s), "", fixed = TRUE)[[1]])
    out <- map[ch]
    out[is.na(out)] <- "N"
    paste0(out, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

#' Convert a sequence to space-separated overlapping k-mers
#'
#' DNABERT consumes sequences as whitespace-delimited k-mers. `k = 6` matches
#' the fine-tuned models in the DeepVRegulome collection.
#'
#' @param seq Character vector of sequences.
#' @param k K-mer size.
#' @return Character vector; `NA` where the sequence is shorter than `k`.
#' @export
#' @examples
#' seq2kmer("ACGTACGT", k = 6)
seq2kmer <- function(seq, k = 6L) {
  k <- as.integer(k)
  vapply(seq, function(s) {
    if (is.na(s) || nchar(s) < k) return(NA_character_)
    n <- nchar(s)
    paste(substring(s, seq_len(n - k + 1L), seq_len(n - k + 1L) + k - 1L),
          collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}
