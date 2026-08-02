## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)
library(PGViS)

## ----stages, eval = TRUE------------------------------------------------------
pgvis_stages()[, c("stage", "description")]

## ----config-------------------------------------------------------------------
#  cfg <- pgvis_config(
#    root            = "/scratch/me/lung_pgvis",
#    reference_fasta = "/refs/hg38.fa",
#    dbsnp_dir       = "/refs/dbsnp155/by_chr",
#    regions   = list(acceptor = "/refs/acceptor_sorted.bed",
#                     donor    = "/refs/donor_sorted.bed"),
#    tfbs_dir  = "/refs/encode_tfbs",
#    cohorts   = c("luad", "lusc"),
#    threads   = 8,
#    scheduler = "slurm",
#    slurm     = list(partition = "cpu", mem = "48G", time = "2-00:00:00",
#                     cpus_per_task = 8, account = "mylab")
#  )
#  pgvis_init(cfg)

## ----sheet--------------------------------------------------------------------
#  sheet <- data.frame(
#    cohort = c("luad", "luad", "lusc"),
#    tissue = c("blood", "solid", "blood"),
#    sample = c("PT-001", "PT-001", "PT-002")
#  )
#  pgvis_sample_sheet(sheet,
#                     path    = file.path(cfg$root, "samples.tsv"),
#                     bam_dir = "/data/raw_bams")
#  #> Warning: 1/3 samples have no BAM under /data/raw_bams. First few: PT-002

## ----run----------------------------------------------------------------------
#  pgvis_merge_bams(cfg, "samples.tsv", submit = TRUE)
#  pgvis_call_gvcf(cfg,  "samples.tsv", submit = TRUE)   # array over chromosomes
#  pgvis_combine_gvcf(cfg,    submit = TRUE)              # array over cohort x chr
#  pgvis_genotype_filter(cfg, submit = TRUE)
#  pgvis_vcf_to_bed(cfg,      submit = TRUE)
#  pgvis_intersect(cfg,       submit = TRUE)              # array over features

## ----inspect------------------------------------------------------------------
#  script <- pgvis_call_gvcf(cfg, "samples.tsv")
#  cat(readLines(script)[1:25], sep = "\n")

## ----force--------------------------------------------------------------------
#  pgvis_intersect(cfg, submit = TRUE, PGVIS_FORCE = 1)

## ----override-----------------------------------------------------------------
#  pgvis_genotype_filter(
#    cfg, submit = TRUE,
#    FILTER_EXPR = "QD < 2.0 || FS > 60.0 || MQ < 40.0",
#    JAVA_MEM    = "32g"
#  )

## ----load---------------------------------------------------------------------
#  variants <- pgvis_load_variants(cfg)
#  pgvis_qc(variants)
#  pgvis_summarise_by(variants, by = c("cohort", "feature_type"))

## ----panel--------------------------------------------------------------------
#  panel <- pgvis_panel_config(
#    root    = "/scratch/me/1kgp_pgvis",
#    vcf_dir = "/refs/1kgp/vcf",
#    regions = cfg$regions,
#    dbsnp_dir = cfg$dbsnp_dir
#  )
#  pgvis_init(panel)
#  
#  pgvis_vcf_to_bed(panel, submit = TRUE)
#  pgvis_intersect(panel,  submit = TRUE)
#  
#  panel_variants <- pgvis_load_variants(panel, region_first = FALSE)
#  ped <- pgvis_load_pedigree("/refs/1kgp/20130606_g1k_3202_samples_ped_population.txt")

## ----bind---------------------------------------------------------------------
#  stacked <- pgvis_bind_panel(cohort_per_sample, panel_per_sample)
#  scores  <- pgvis_score(stacked, reference_group = "1000Genomes")

## ----manifest-----------------------------------------------------------------
#  pgvis_manifest(cfg, stage = "custom_filter",
#                 inputs  = list(n_in = nrow(variants)),
#                 outputs = list(rds = "filtered.rds"))

