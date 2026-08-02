## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PGViS)
set.seed(1)

## ----metrics-demo-------------------------------------------------------------
# strongly disruptive: model is confident before, not after
log_odds_ratio(p_ref = 0.97, p_alt = 0.03)
score_change(0.97, 0.03)

# a variant enriched in the cohort relative to the population
d <- data.frame(af_ref = 0.60, af_alt = 0.40,
                pop_af_ref = 0.95, pop_af_alt = 0.05)
pgvis_delta_wi(d, transform = FALSE, keep_raw = TRUE)[, c("delta_af_alt", "delta_wi")]

## ----shrinkage----------------------------------------------------------------
c(confident = score_change(0.9, 0.7), unconfident = score_change(0.3, 0.1))

## ----config-------------------------------------------------------------------
cfg <- pgvis_config(
  root       = file.path(tempdir(), "demo"),
  cohorts    = c("luad", "lusc"),
  thresholds = list(ref_prob_min = 0.7, alt_prob_max = 0.3)
)
cfg

## ----stages-------------------------------------------------------------------
pgvis_stages()[, c("stage", "description")]

## ----render, eval = FALSE-----------------------------------------------------
#  pgvis_call_gvcf(cfg, sample_sheet = "samples.tsv", submit = FALSE)
#  #> [PGViS] Rendered 02_call_gvcf.sh -> .../logs/02_call_gvcf.run.sh

## ----keys---------------------------------------------------------------------
k <- func_site_key("chr7", 55019000, 55019080, "-")
v <- variant_key("chr7", 55019032, 55019033, "a", "g")
c(k, v)
parse_key(v, "variant")

## ----keys-strandless----------------------------------------------------------
func_site_key("chr7", 100, 400, ".", stranded = FALSE)

## ----mutate-------------------------------------------------------------------
mutate_sequence("A", "G", "TTTATTT", pos_start = 4, pos_end = 4)  # exact
mutate_sequence("A", "G", "TTTATTT", pos_start = 5, pos_end = 5)  # off by one
mutate_sequence("CCCC", "G", "TTTATTT", pos_start = 4, pos_end = 7)  # wrong allele

## ----kmer---------------------------------------------------------------------
seq2kmer("ACGTACGTAC", k = 6)

## ----popaf--------------------------------------------------------------------
row <- data.frame(
  dbsnp_freqs = paste(c("0.11", "", "0.13", rep("", 9), "0.19"), collapse = ","),
  dbsnp_minor = paste(c("A",    "", "A",    rep("", 9), "A"),    collapse = ","),
  dbsnp_ref = "G", dbsnp_alt = "A"
)
pgvis_resolve_population_af(row)[, c("pop_af_alt", "pop_af_ref", "freq_source")]

## ----popaf-priority-----------------------------------------------------------
pgvis_resolve_population_af(row, priority = c("1000G", "TOPMED"))$freq_source

## ----score--------------------------------------------------------------------
n <- 150
per_sample <- rbind(
  data.frame(sample = paste0("R", 1:n), cohort = "1000Genomes",
             superpopulation = rep(c("EUR", "EAS"), length.out = n),
             n_variants = 10L,
             delta_wi_int = rnorm(n), log_odds_ratio = rnorm(n),
             delta_af_alt = rnorm(n)),
  data.frame(sample = paste0("T", 1:n), cohort = "tumour",
             superpopulation = rep(c("EUR", "EAS"), length.out = n),
             n_variants = 13L,
             delta_wi_int = rnorm(n, 1.2), log_odds_ratio = rnorm(n, 1.2),
             delta_af_alt = rnorm(n, 1.2))
)

scored <- pgvis_score(per_sample, reference_group = "1000Genomes")
head(scored[, c("sample", "cohort", "superpopulation", "pgvis")])
attr(scored, "weights")

## ----score-check--------------------------------------------------------------
tapply(scored$pgvis, scored$cohort, median)

## ----frozen, eval = FALSE-----------------------------------------------------
#  w  <- attr(discovery_scores, "weights")
#  val <- pgvis_score(validation, weights = w)

## ----compare------------------------------------------------------------------
d <- data.frame(pgvis = c(rnorm(120, 0), rnorm(120, 1.1)),
                group = rep(c("control", "case"), each = 120))
pgvis_compare_groups(d, min_n = 50, n_perm = 500)[, c("n1", "n2", "auc", "cohens_d", "perm_p")]

## ----compare-small------------------------------------------------------------
small <- data.frame(pgvis = rnorm(20), group = rep(c("a", "b"), 10))
pgvis_compare_groups(small, min_n = 50)$label

## ----si-----------------------------------------------------------------------
sessionInfo()

