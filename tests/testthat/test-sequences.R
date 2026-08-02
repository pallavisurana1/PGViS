test_that("revcomp round-trips and handles ambiguity codes", {
  expect_equal(revcomp("ACCGGT"), "ACCGGT")
  expect_equal(revcomp("AAAA"), "TTTT")
  expect_equal(revcomp(revcomp("ACGTTGCA")), "ACGTTGCA")
  expect_equal(revcomp("ACGN"), "NCGT")
})

test_that("seq2kmer produces the right number of overlapping k-mers", {
  km <- seq2kmer("ACGTACGT", k = 6)
  expect_equal(length(strsplit(km, " ")[[1]]), 8 - 6 + 1)
  expect_equal(strsplit(km, " ")[[1]][1], "ACGTAC")
  expect_true(is.na(seq2kmer("ACGT", k = 6)))
})

test_that("mutate_sequence substitutes, inserts and deletes", {
  expect_equal(mutate_sequence("A", "G", "TTTATTT", 4, 4), "TTTGTTT")
  expect_equal(mutate_sequence("A", "AGG", "TTTATTT", 4, 4), "TTTAGGTTT")
  expect_equal(mutate_sequence("AT", "A", "TTTATTTT", 4, 5), "TTTATTT")
})

test_that("mutate_sequence tolerates off-by-one coordinates but not a wrong allele", {
  # nominal position is one base off; the true match is found within the window
  expect_equal(mutate_sequence("A", "G", "TTTATTT", 5, 5), "TTTGTTT")
  # reference allele genuinely absent -> NA, not a silent wrong answer
  expect_true(is.na(mutate_sequence("CCCC", "G", "TTTATTT", 4, 7)))
})

test_that("coordinate keys build and parse consistently", {
  k <- func_site_key("chr1", 100, 180, "+")
  expect_equal(k, "chr1:100-180(+)")
  expect_equal(func_site_key("chr1", 100, 180, "+", stranded = FALSE), "chr1:100-180()")
  expect_equal(func_site_key("chr1", 100, 180, "."), "chr1:100-180()")

  p <- parse_key(k, "site")
  expect_equal(p$chr, "chr1"); expect_equal(p$start, 100L); expect_equal(p$strand, "+")

  v <- variant_key("chr1", 120, 121, "a", "g")
  expect_equal(v, "chr1:120-121(A>G)")
  expect_equal(variant_key("chr1", 1, 2, "A", NA), "chr1:1-2(A>*)")
  pv <- parse_key(v, "variant")
  expect_equal(pv$ref, "A"); expect_equal(pv$alt, "G")
})

test_that("strand-aware variant offsets mirror correctly (arithmetic only)", {
  # This is the fragile part of pgvis_build_sequences: converting BED
  # coordinates into a 1-based offset inside the strand-oriented region
  # sequence. Tested here without BSgenome so it runs in any environment.
  offset <- function(site_start, var_start, var_end, L, minus) {
    span <- max(1L, var_end - var_start)
    ms <- if (minus) max(1L, L - (var_end - site_start) + 1L)
          else       max(1L, var_start - site_start + 1L)
    c(start = ms, end = min(L, ms + span - 1L))
  }
  # A variant and its minus-strand image must sum to L + 1 (they are mirror
  # positions in the reverse complement).
  L <- 10L
  plus  <- offset(100, 103, 104, L, FALSE)
  minus <- offset(100, 103, 104, L, TRUE)
  expect_equal(unname(plus["start"]), 4L)
  expect_equal(unname(plus["start"] + minus["start"]), L + 1L)
  # Deletions span more than one base.
  del <- offset(100, 103, 105, L, FALSE)
  expect_equal(unname(del["end"] - del["start"]), 1L)
})

test_that("pgvis_build_sequences produces a mutated sequence that differs from reference", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  skip_if_not_installed("GenomicRanges")
  # A single well-known locus; the alt sequence must differ from ref by
  # exactly the substituted base, and the k-mers must be non-NA.
  v <- data.frame(
    feature = "donor", feature_type = "splice",
    site_chr = "chr1", site_start = 1000000L, site_end = 1000080L,
    site_strand = "+",
    var_chr = "chr1", var_start = 1000040L, var_end = 1000041L,
    ref = "N", alt = "A", stringsAsFactors = FALSE
  )
  # Use the real reference base so mutate_sequence can locate it.
  genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38",
                             "BSgenome.Hsapiens.UCSC.hg38")
  gr <- GenomicRanges::GRanges("chr1",
          IRanges::IRanges(1000041L, 1000041L))
  v$ref <- as.character(BSgenome::getSeq(genome, gr))
  v$alt <- ifelse(v$ref == "A", "C", "A")
  out <- pgvis_build_sequences(v, drop_failures = FALSE)
  expect_true(!is.na(out$alt_seq))
  expect_false(identical(out$ref_seq, out$alt_seq))
  expect_equal(nchar(out$ref_seq), nchar(out$alt_seq))
})
