## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PGViS)

## ----lor----------------------------------------------------------------------
log_odds_ratio(p_ref = 0.97, p_alt = 0.03)
log_odds_ratio(0.5, 0.5)          # no change -> 0

## ----sc-----------------------------------------------------------------------
score_change(0.9, 0.7)    # confident, so the shift counts
score_change(0.3, 0.1)    # same 0.2 shift, but low confidence, so shrunk

## ----dwi----------------------------------------------------------------------
d <- data.frame(af_ref = 0.60, af_alt = 0.40,
                pop_af_ref = 0.95, pop_af_alt = 0.05)
pgvis_delta_wi(d, transform = FALSE, keep_raw = TRUE)[, c("delta_af_alt", "delta_wi")]

## ----popaf--------------------------------------------------------------------
row <- data.frame(
  dbsnp_freqs = paste(c("0.11", "", "0.13", rep("", 9), "0.19"), collapse = ","),
  dbsnp_minor = paste(c("A",    "", "A",    rep("", 9), "A"),    collapse = ","),
  dbsnp_ref = "G", dbsnp_alt = "A"
)
pgvis_resolve_population_af(row)[, c("pop_af_alt", "pop_af_ref", "freq_source")]

## ----popaf-priority-----------------------------------------------------------
pgvis_resolve_population_af(row, priority = c("1000G", "TOPMED"))$freq_source

## ----metrics------------------------------------------------------------------
variants <- pgvis_example("variants")
preds    <- pgvis_example("predictions")

metrics <- pgvis_variant_metrics(variants, preds)
metrics[1:3, c("feature", "log_odds_ratio", "score_change", "delta_wi_int")]

## ----filter-------------------------------------------------------------------
kept <- pgvis_filter_disruptive(
  metrics,
  ref_prob_min = 0.5, alt_prob_max = 0.5, min_carrier_frac = 0
)
nrow(metrics)
nrow(kept)

## ----audit--------------------------------------------------------------------
attr(kept, "pgvis_filter")

## ----qc-----------------------------------------------------------------------
pgvis_qc(metrics)

## ----replication--------------------------------------------------------------
combined <- data.frame(
  study     = c("tcga", "cptac", "tcga"),
  cancer    = "luad",
  func_site = c("f1", "f1", "f2"),
  var_site  = c("v1", "v1", "v2"),
  stringsAsFactors = FALSE
)
pgvis_filter_replicated(combined, within = "cancer", min_studies = 2)

