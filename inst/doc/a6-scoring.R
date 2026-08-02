## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PGViS)
set.seed(1)

## ----explode------------------------------------------------------------------
metrics <- pgvis_variant_metrics(pgvis_example("variants"),
                                 pgvis_example("predictions"))
metrics <- pgvis_filter_disruptive(metrics, ref_prob_min = 0.5,
                                   alt_prob_max = 0.5, min_carrier_frac = 0)

long <- pgvis_explode_carriers(metrics)
nrow(metrics)
nrow(long)          # one row per (variant, carrier)

## ----aggregate----------------------------------------------------------------
clin <- pgvis_example("clinical")
clin$superpopulation <- pgvis_map_superpopulation(clin$race)
long <- dplyr::left_join(long, clin[, c("sample", "superpopulation")],
                         by = "sample")

per_sample <- pgvis_aggregate_samples(long)
per_sample <- per_sample[!is.na(per_sample$superpopulation), ]
head(per_sample)

## ----panel--------------------------------------------------------------------
ped <- pgvis_load_pedigree(
  system.file("extdata", "toy_pedigree.txt", package = "PGViS"))

panel <- data.frame(
  sample = ped$sample, superpopulation = ped$superpopulation,
  cohort = "1000Genomes", n_variants = 10L,
  delta_wi_int   = rnorm(nrow(ped)),
  log_odds_ratio = rnorm(nrow(ped)),
  delta_af_alt   = rnorm(nrow(ped))
)

stacked <- pgvis_bind_panel(per_sample, panel)

## ----score--------------------------------------------------------------------
scored <- pgvis_score(stacked, reference_group = "1000Genomes")
head(scored[, c("sample", "cohort", "superpopulation", "pgvis")])

## ----weights------------------------------------------------------------------
attr(scored, "weights")

## ----frozen-------------------------------------------------------------------
w  <- attr(scored, "weights")
s2 <- pgvis_score(stacked, weights = w)      # same scale as `scored`
identical(attr(s2, "weights"), w)

## ----compare------------------------------------------------------------------
demo <- data.frame(pgvis = c(rnorm(120, 0), rnorm(120, 1.1)),
                   group = rep(c("control", "case"), each = 120))
pgvis_compare_groups(demo, min_n = 50, n_perm = 500)[,
  c("n1", "n2", "auc", "cohens_d", "perm_p")]

## ----survival-----------------------------------------------------------------
n <- 300
score <- rnorm(n)
time  <- rexp(n, rate = exp(0.7 * score) / 500)   # higher score, shorter time
sc  <- data.frame(sample = paste0("S", 1:n), pgvis = score)
cl  <- data.frame(sample = paste0("S", 1:n), time = time,
                  event = rbinom(n, 1, 0.8))

res <- pgvis_survival(sc, cl)
res$summary

## ----run, eval = FALSE--------------------------------------------------------
#  res <- pgvis_run(cfg, variants, clinical = clin,
#                   reference_group = "1000Genomes")

