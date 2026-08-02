#' Create a PGViS project configuration
#'
#' Every path, threshold and external-tool location used anywhere in the
#' workflow lives in a single `pgvis_config` object. This replaces the
#' hard-coded absolute paths (`/vast/home/...`) and per-script
#' machine-specific settings, and it is what makes a project portable to
#' another machine without editing code.
#'
#' The directory layout under `root` is created on demand and is stable:
#'
#' ```
#' root/
#'   01_bam/          input / merged BAMs
#'   02_gvcf/         per-sample per-chromosome GVCFs
#'   03_joint/        joint-genotyped, filtered, normalised cohort VCFs
#'   04_bed/          cohort BED + carrier annotation
#'   05_intersect/    regulatory-region intersections
#'   06_dbsnp/        dbSNP-annotated intersections
#'   07_seq/          ref/alt sequences and DNABERT k-mer inputs
#'   08_preds/        model predictions
#'   09_metrics/      per-variant metric tables
#'   10_score/        per-sample PGViS scores and figures
#'   logs/
#' ```
#'
#' @param root Project root directory. Created if absent.
#' @param genome Genome build label. Only `"hg38"` is currently supported.
#' @param reference_fasta Path to the reference FASTA used for variant calling.
#' @param dbsnp_dir Directory of per-chromosome dbSNP BED files
#'   (e.g. `dbSnp155.chr1.bed`). See [pgvis_prepare_dbsnp()].
#' @param regions A named list of BED files defining the regulatory regions to
#'   test. Names become the `feature` column downstream. Typically
#'   `list(acceptor = "...", donor = "...")` plus one entry per TFBS model.
#' @param tfbs_dir Optional directory containing one subdirectory per TF, each
#'   holding a region BED. Expanded into `regions` when supplied.
#' @param chromosomes Chromosomes to process.
#' @param cohorts Character vector of cohort labels (e.g. `c("luad","lusc")`).
#'   Replaces the duplicated `*_LUAD.sh` / `*_LUSC.sh` scripts.
#' @param threads Threads for external tools.
#' @param scheduler One of `"none"`, `"slurm"`. Controls how
#'   [pgvis_submit()] launches the shell stages.
#' @param slurm A named list of `#SBATCH` defaults (partition, mem, time, ...).
#' @param python Optional path to a Python executable / conda env that has
#'   `deepvregulome` installed. `NULL` means "resolve at call time".
#' @param thresholds Named list of analysis thresholds. See Details.
#'
#' @details
#' `thresholds` controls the filtering applied downstream:
#' \describe{
#'   \item{`ref_prob_min`}{Minimum model probability on the reference
#'     sequence for a site to be considered a genuine functional element
#'     (0.5 in the permissive setting, 0.7 in the strict one).}
#'   \item{`alt_prob_max`}{Maximum probability on the mutated sequence for the
#'     variant to be called disruptive (0.5 permissive, 0.3 strict).}
#'   \item{`min_carrier_frac`}{Minimum carriers as a fraction of the maximum
#'     observed carrier count (0.05 = the "5 percent of patients" filter).}
#'   \item{`max_indel_len`}{Longest REF/ALT allele kept (10).}
#'   \item{`pseudocount`}{Added to all frequencies before ratios (1e-6).}
#'   \item{`drop_freq_sources`}{Population-frequency sources to exclude when
#'     scoring, so that a 1000G reference panel is not compared against
#'     1000G-derived frequencies.}
#' }
#'
#' @return An object of class `pgvis_config`.
#' @export
#' @examples
#' cfg <- pgvis_config(root = tempfile(), reference_fasta = NA_character_)
#' cfg$dirs$score
pgvis_config <- function(root,
                         genome = "hg38",
                         reference_fasta = NULL,
                         dbsnp_dir = NULL,
                         regions = list(),
                         tfbs_dir = NULL,
                         chromosomes = c(as.character(1:22), "X", "Y"),
                         cohorts = character(),
                         threads = 4L,
                         scheduler = c("none", "slurm"),
                         slurm = list(partition = NULL, mem = "32G",
                                      time = "24:00:00", cpus_per_task = 4L,
                                      account = NULL),
                         python = NULL,
                         thresholds = list()) {

  scheduler <- match.arg(scheduler)
  stopifnot(is.character(root), length(root) == 1L)
  if (!identical(genome, "hg38")) {
    stop("Only genome = 'hg38' is currently supported.", call. = FALSE)
  }

  default_thresholds <- list(
    ref_prob_min     = 0.5,
    alt_prob_max     = 0.5,
    min_carrier_frac = 0.05,
    max_indel_len    = 10L,
    pseudocount      = 1e-6,
    drop_freq_sources = "1000G"
  )
  thresholds <- utils::modifyList(default_thresholds, thresholds)

  if (!is.null(tfbs_dir)) {
    regions <- c(regions, .pgvis_expand_tfbs_dir(tfbs_dir))
  }
  if (length(regions) && is.null(names(regions))) {
    stop("`regions` must be a *named* list of BED files.", call. = FALSE)
  }

  subdirs <- c(bam = "01_bam", gvcf = "02_gvcf", joint = "03_joint",
               bed = "04_bed", intersect = "05_intersect", dbsnp = "06_dbsnp",
               seq = "07_seq", preds = "08_preds", metrics = "09_metrics",
               score = "10_score", logs = "logs")
  dirs <- as.list(file.path(root, subdirs))
  names(dirs) <- names(subdirs)

  structure(
    list(root = root, genome = genome, dirs = dirs,
         reference_fasta = reference_fasta, dbsnp_dir = dbsnp_dir,
         regions = regions, chromosomes = chromosomes, cohorts = cohorts,
         threads = as.integer(threads), scheduler = scheduler, slurm = slurm,
         python = python, thresholds = thresholds),
    class = "pgvis_config"
  )
}

#' @export
print.pgvis_config <- function(x, ...) {
  cat("<pgvis_config>\n")
  cat("  root       :", x$root, "\n")
  cat("  genome     :", x$genome, "\n")
  cat("  cohorts    :", if (length(x$cohorts)) paste(x$cohorts, collapse = ", ") else "-", "\n")
  cat("  chromosomes:", length(x$chromosomes), "\n")
  cat("  regions    :", length(x$regions), "feature BED file(s)\n")
  cat("  scheduler  :", x$scheduler, "\n")
  cat("  thresholds : ref_prob >", x$thresholds$ref_prob_min,
      "| alt_prob <=", x$thresholds$alt_prob_max,
      "| carriers >=", paste0(100 * x$thresholds$min_carrier_frac, "%"), "\n")
  invisible(x)
}

#' Create the on-disk directory tree for a project
#'
#' @param cfg A [pgvis_config()].
#' @return `cfg`, invisibly.
#' @export
pgvis_init <- function(cfg) {
  .assert_config(cfg)
  for (d in unlist(cfg$dirs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(cfg)
}

.pgvis_expand_tfbs_dir <- function(tfbs_dir, pattern = "\\.bed$") {
  subs <- list.dirs(tfbs_dir, recursive = FALSE, full.names = TRUE)
  out <- lapply(subs, function(d) {
    f <- list.files(d, pattern = pattern, full.names = TRUE)
    if (length(f)) f[[1]] else NULL
  })
  names(out) <- basename(subs)
  out[!vapply(out, is.null, logical(1))]
}

.assert_config <- function(cfg) {
  if (!inherits(cfg, "pgvis_config")) {
    stop("Expected a `pgvis_config` object; see ?pgvis_config.", call. = FALSE)
  }
  invisible(TRUE)
}
