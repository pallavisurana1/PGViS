# ---------------------------------------------------------------------------
# Reference data.
#
# The pipeline needs three external reference resources: the genome FASTA, the
# dbSnp155 frequency track, and the regulatory-region BEDs. This module
# documents where each comes from and can fetch and prepare the dbSNP track,
# which is the one that needs real work: it ships as a single binary bigBed and
# has to be split per chromosome to match how the intersection stage reads it.
# ---------------------------------------------------------------------------

#' Canonical URLs for PGViS reference data
#'
#' All point at the UCSC Genome Browser download server, which is the source
#' the dbSnp155 frequency parsing in [pgvis_resolve_population_af()] was built
#' against. Using a different dbSNP build, or the NCBI VCF rather than the UCSC
#' track, will shift the comma-delimited frequency field and silently break
#' the source-index lookup -- see [pgvis_freq_sources()].
#'
#' @param genome Genome build. Only `"hg38"` is supported.
#' @return A named list of URLs.
#' @export
#' @examples
#' pgvis_reference_urls()$dbsnp155_all
pgvis_reference_urls <- function(genome = "hg38") {
  if (!identical(genome, "hg38")) {
    stop("Only hg38 reference data is wired up.", call. = FALSE)
  }
  base <- "https://hgdownload.soe.ucsc.edu"
  list(
    # dbSnp155: full track (~1 billion variants) and the common subset.
    # Both are bigBed and must be converted with bigBedToBed.
    dbsnp155_all    = paste0(base, "/gbdb/hg38/snp/dbSnp155.bb"),
    dbsnp155_common = paste0(base, "/gbdb/hg38/snp/dbSnp155Common.bb"),
    # UCSC bigBed utilities (Linux x86_64). macOS: swap linux.x86_64 for macOSX.x86_64.
    bigbedtobed     = paste0(base, "/admin/exe/linux.x86_64/bigBedToBed"),
    bigbedinfo      = paste0(base, "/admin/exe/linux.x86_64/bigBedInfo"),
    # NCBI analysis-set FASTA: the build recommended for variant calling
    # (alt/fix sequences removed, PARs hard-masked).
    genome_fasta    = paste0(base, "/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"),
    chrom_sizes     = paste0(base, "/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chrom.sizes")
  )
}

#' Prepare the per-chromosome dbSNP frequency track
#'
#' Downloads the dbSnp155 bigBed, extracts it to BED with UCSC's `bigBedToBed`,
#' and splits it into the per-chromosome files (`dbSnp155.chr1.bed`, ...) that
#' the intersection stage expects in `cfg$dbsnp_dir`.
#'
#' This is a large download (the full track is tens of GB extracted). For most
#' analyses the common subset (`subset = "common"`) is enough and far smaller.
#' The function is restartable: chromosomes whose output already exists and is
#' non-empty are skipped.
#'
#' @param dest Output directory. Set `cfg$dbsnp_dir` to this.
#' @param subset `"common"` (recommended, ~a few GB) or `"all"`.
#' @param chromosomes Which chromosomes to produce.
#' @param bigbedtobed Path to the `bigBedToBed` binary. If `NULL`, PGViS looks
#'   on `PATH` and, failing that, tells you where to get it.
#' @param genome Genome build.
#' @param keep_bigbed Keep the downloaded `.bb` afterwards.
#' @param dry_run Print the commands instead of running them. Useful on a
#'   login node where you want to inspect before launching a long job.
#' @return The output directory, invisibly.
#' @export
#' @examples
#' # Show what it would do without downloading anything:
#' pgvis_prepare_dbsnp(tempdir(), subset = "common", dry_run = TRUE)
pgvis_prepare_dbsnp <- function(dest,
                                subset = c("common", "all"),
                                chromosomes = c(as.character(1:22), "X", "Y"),
                                bigbedtobed = NULL,
                                genome = "hg38",
                                keep_bigbed = FALSE,
                                dry_run = FALSE) {
  subset <- match.arg(subset)
  urls <- pgvis_reference_urls(genome)
  bb_url <- if (subset == "common") urls$dbsnp155_common else urls$dbsnp155_all
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  bb_local <- file.path(dest, basename(bb_url))
  chroms <- paste0(ifelse(grepl("^chr", chromosomes), "", "chr"), chromosomes)

  # Resolve the bigBedToBed binary.
  tool <- bigbedtobed %||% unname(Sys.which("bigBedToBed"))
  if (!nzchar(tool) && !dry_run) {
    stop("`bigBedToBed` not found. Download the UCSC utility from\n    ",
         urls$bigbedtobed, "\n",
         "make it executable, and pass its path as `bigbedtobed`, ",
         "or add it to PATH.", call. = FALSE)
  }

  steps <- c(
    sprintf("# 1. Download the dbSnp155 %s track (%s)", subset,
            if (subset == "common") "smaller, recommended" else "very large"),
    sprintf("curl -o %s %s", shQuote(bb_local), shQuote(bb_url)),
    "# 2. Extract each chromosome directly from the bigBed (no full extraction needed)"
  )
  tool_display <- if (nzchar(tool %||% "")) tool else "bigBedToBed"
  for (chr in chroms) {
    out <- file.path(dest, sprintf("dbSnp155.%s.bed", chr))
    steps <- c(steps, sprintf("%s -chrom=%s %s %s",
                              shQuote(tool_display), chr,
                              shQuote(bb_local), shQuote(out)))
  }
  if (!keep_bigbed) steps <- c(steps, sprintf("rm -f %s", shQuote(bb_local)))

  if (dry_run) {
    cat(steps, sep = "\n"); cat("\n")
    return(invisible(dest))
  }

  .msg("Downloading ", bb_url)
  if (!file.exists(bb_local) || file.size(bb_local) == 0) {
    ok <- utils::download.file(bb_url, bb_local, mode = "wb", quiet = FALSE)
    if (ok != 0) stop("Download failed. Check network access to hgdownload.soe.ucsc.edu.",
                      call. = FALSE)
  }

  # bigBedToBed can subset by chromosome without extracting the whole file,
  # which is the difference between minutes and hours on the full track.
  for (chr in chroms) {
    out <- file.path(dest, sprintf("dbSnp155.%s.bed", chr))
    if (file.exists(out) && file.size(out) > 0) {
      .msg("skip ", chr, " (exists)"); next
    }
    .msg("extracting ", chr)
    st <- system2(tool, c(sprintf("-chrom=%s", chr),
                          shQuote(bb_local), shQuote(out)))
    if (st != 0) warning("bigBedToBed failed for ", chr, call. = FALSE)
  }
  if (!keep_bigbed) unlink(bb_local)

  n <- length(list.files(dest, pattern = "^dbSnp155\\.chr.*\\.bed$"))
  .msg(sprintf("dbSNP prep done: %d per-chromosome BED file(s) in %s", n, dest))
  invisible(dest)
}

#' Check that a project's reference data is present and usable
#'
#' Run this before launching the pipeline. It verifies the reference FASTA is
#' indexed, the dbSNP directory has per-chromosome files for the configured
#' chromosomes, and the region BEDs exist and are non-empty -- the three
#' things whose absence otherwise surfaces as a cryptic failure several hours
#' into a cluster job.
#'
#' @param cfg A [pgvis_config()].
#' @return A data frame of checks with a `pass` column, invisibly. Prints a
#'   readable report.
#' @export
pgvis_check_references <- function(cfg) {
  .assert_config(cfg)
  checks <- list()
  add <- function(what, pass, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = what, pass = pass, detail = detail, stringsAsFactors = FALSE)
  }

  # Reference FASTA + index.
  fa <- cfg$reference_fasta
  if (is.null(fa) || is.na(fa) || !nzchar(fa)) {
    add("reference FASTA", FALSE, "cfg$reference_fasta is not set")
  } else if (!file.exists(fa)) {
    add("reference FASTA", FALSE, paste("missing:", fa))
  } else {
    add("reference FASTA", TRUE, fa)
    fai <- paste0(fa, ".fai")
    add("FASTA index (.fai)", file.exists(fai),
        if (file.exists(fai)) fai else "run: samtools faidx <fasta>")
  }

  # dbSNP per-chromosome files.
  if (is.null(cfg$dbsnp_dir) || !nzchar(cfg$dbsnp_dir)) {
    add("dbSNP directory", FALSE, "cfg$dbsnp_dir is not set (see pgvis_prepare_dbsnp)")
  } else {
    want <- paste0("dbSnp155.chr",
                   sub("^chr", "", cfg$chromosomes), ".bed")
    have <- file.exists(file.path(cfg$dbsnp_dir, want))
    add("dbSNP per-chromosome BEDs", all(have),
        sprintf("%d/%d present in %s", sum(have), length(have), cfg$dbsnp_dir))
  }

  # Region BEDs.
  if (!length(cfg$regions)) {
    add("region BEDs", FALSE, "cfg$regions is empty")
  } else {
    exist <- vapply(cfg$regions, function(f)
      !is.null(f) && file.exists(f) && file.size(f) > 0, logical(1))
    add("region BEDs", all(exist),
        sprintf("%d/%d present and non-empty", sum(exist), length(exist)))
  }

  res <- do.call(rbind, checks)
  cat("PGViS reference check\n---------------------\n")
  for (i in seq_len(nrow(res))) {
    cat(sprintf("[%s] %-28s %s\n",
                if (res$pass[i]) "OK" else "!!", res$check[i], res$detail[i]))
  }
  if (!all(res$pass)) {
    cat("\nSome checks failed; the pipeline will not run cleanly until they pass.\n")
  }
  invisible(res)
}
