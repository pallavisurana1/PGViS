## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PGViS)

## ----config-------------------------------------------------------------------
cfg <- pgvis_config(
  root      = file.path(tempdir(), "demo_project"),
  cohorts   = c("luad", "lusc")
)
cfg

## ----dirs---------------------------------------------------------------------
cfg$dirs$gvcf
cfg$dirs$score

## ----init---------------------------------------------------------------------
pgvis_init(cfg)
dir.exists(cfg$dirs$score)

## ----thresholds---------------------------------------------------------------
strict <- pgvis_config(
  root = tempdir(),
  thresholds = list(ref_prob_min = 0.7, alt_prob_max = 0.3)
)
strict$thresholds

## ----urls---------------------------------------------------------------------
urls <- pgvis_reference_urls("hg38")
urls$dbsnp155_common
urls$genome_fasta

## ----dbsnp-dry----------------------------------------------------------------
pgvis_prepare_dbsnp(
  dest        = file.path(cfg$root, "dbsnp"),
  subset      = "common",          # far smaller than "all"; enough for most work
  chromosomes = c("21", "22"),
  dry_run     = TRUE
)

## ----check--------------------------------------------------------------------
pgvis_check_references(cfg)

