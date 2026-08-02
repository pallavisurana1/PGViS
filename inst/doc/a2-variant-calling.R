## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PGViS)

## ----stages-------------------------------------------------------------------
pgvis_stages()[, c("stage", "description")]

## ----config, eval = FALSE-----------------------------------------------------
#  cfg <- pgvis_config(
#    root            = "/scratch/me/lung_pgvis",
#    reference_fasta = "/refs/hg38.analysisSet.fa",
#    dbsnp_dir       = "/refs/dbsnp155/by_chr",
#    regions   = list(acceptor = "/refs/acceptor.bed",
#                     donor    = "/refs/donor.bed"),
#    tfbs_dir  = "/refs/encode_tfbs",          # one subdirectory per TF
#    cohorts   = c("luad", "lusc"),
#    threads   = 8,
#    scheduler = "slurm",
#    slurm     = list(partition = "cpu", mem = "48G", time = "2-00:00:00",
#                     cpus_per_task = 8, account = "mylab")
#  )
#  pgvis_init(cfg)

## ----sheet--------------------------------------------------------------------
sheet <- data.frame(
  cohort = c("luad", "luad", "lusc"),
  tissue = c("blood", "solid", "blood"),
  sample = c("PT-001", "PT-001", "PT-002")
)
pgvis_sample_sheet(sheet)

## ----run, eval = FALSE--------------------------------------------------------
#  pgvis_merge_bams(cfg, "samples.tsv", submit = TRUE)   # 01: merge replicates
#  pgvis_call_gvcf(cfg,  "samples.tsv", submit = TRUE)   # 02: HaplotypeCaller, array over chr
#  pgvis_combine_gvcf(cfg,    submit = TRUE)             # 03: GenomicsDBImport, array over cohort x chr
#  pgvis_genotype_filter(cfg, submit = TRUE)             # 04: genotype, hard-filter, normalise
#  pgvis_vcf_to_bed(cfg,      submit = TRUE)             # 05: BED + carriers
#  pgvis_intersect(cfg,       submit = TRUE)             # 06: intersect regions, then dbSNP

## ----inspect, eval = FALSE----------------------------------------------------
#  script <- pgvis_call_gvcf(cfg, "samples.tsv")
#  cat(readLines(script), sep = "\n")

## ----override, eval = FALSE---------------------------------------------------
#  pgvis_genotype_filter(
#    cfg, submit = TRUE,
#    FILTER_EXPR = "QD < 2.0 || FS > 60.0 || MQ < 40.0",
#    JAVA_MEM    = "32g"
#  )

## ----force, eval = FALSE------------------------------------------------------
#  pgvis_intersect(cfg, submit = TRUE, PGVIS_FORCE = 1)

## ----schema-------------------------------------------------------------------
pgvis_schema("splice", with_dbsnp = TRUE)[1:12]

## ----load, eval = FALSE-------------------------------------------------------
#  variants <- pgvis_load_variants(cfg)
#  pgvis_summarise_by(variants, by = c("cohort", "feature_type"))

## ----fixture------------------------------------------------------------------
variants <- pgvis_example("variants")
str(variants[, c("feature", "site_chr", "var_start", "ref", "alt",
                 "af_alt", "n_carrier", "func_site", "var_site")])

