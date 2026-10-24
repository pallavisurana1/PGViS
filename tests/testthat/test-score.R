set.seed(1)

toy <- function(n = 120) {
  rbind(
    data.frame(sample = paste0("R", 1:n), cohort = "1000Genomes",
               superpopulation = rep(c("EUR", "EAS"), length.out = n),
               n_variants = 10L,
               delta_wi_int   = rnorm(n, 0, 1),
               log_odds_ratio = rnorm(n, 0, 1),
               delta_af_alt   = rnorm(n, 0, 1)),
    data.frame(sample = paste0("C", 1:n), cohort = "tumour",
               superpopulation = rep(c("EUR", "EAS"), length.out = n),
               n_variants = 12L,
               delta_wi_int   = rnorm(n, 1.5, 1),
               log_odds_ratio = rnorm(n, 1.5, 1),
               delta_af_alt   = rnorm(n, 1.5, 1))
  )
}

# --- the percentile / geometric mean / root mean square score -----------------

# Long (sample x variant) table. Each cohort has its own variant list; every
# sample carries between 5 and 30 of its list's variants. With `enrich = TRUE`
# the tumour samples preferentially carry the high-ranking variants of their
# list, which is what the score is designed to detect (percentiles are taken
# within each list, so absolute feature values do not matter, only which
# variants of the list a sample carries).
toy_long <- function(n = 60, enrich = TRUE, n_list_ref = 80, n_list_tum = 200) {
  mk <- function(cohort, n_list, pre, enrich, ids) {
    vars <- data.frame(var_site = paste0(pre, seq_len(n_list)),
                       delta_wi_int = rnorm(n_list),
                       log_odds_ratio = rnorm(n_list),
                       delta_af_alt = rnorm(n_list))
    prob <- if (enrich) rank(abs(vars$delta_wi_int) * rank(vars$log_odds_ratio)) else NULL
    do.call(rbind, lapply(ids, function(s) {
      v <- vars[sample(n_list, sample(5:30, 1), prob = prob), ]
      v$sample <- s; v$cohort <- cohort
      v$superpopulation <- if (as.integer(sub("\\D", "", s)) %% 2) "EUR" else "EAS"
      v
    }))
  }
  rbind(mk("1000Genomes", n_list_ref, "r", FALSE, paste0("R", 1:n)),
        mk("tumour", n_list_tum, "t", enrich, paste0("C", 1:n)))
}

test_that("variant scores are percentiles within each dataset's own list", {
  long <- data.frame(cohort = "a", var_site = c("v1", "v2", "v3", "v4"),
                     delta_wi_int = c(-4, 1, 2, 3), log_odds_ratio = c(1, 2, 3, 4))
  vs <- pgvis_variant_score(long)
  # |dWi| ranks: 4 -> 1, 1 -> .25, 2 -> .5, 3 -> .75 ; LOR ranks .25 .5 .75 1
  expect_equal(vs$p_delta_wi_int, c(1, .25, .5, .75))
  expect_equal(vs$p_log_odds_ratio, c(.25, .5, .75, 1))
  expect_equal(vs$variant_score, sqrt(vs$p_delta_wi_int * vs$p_log_odds_ratio))
  expect_true(all(vs$variant_score > 0 & vs$variant_score <= 1))
})

test_that("percentiles are computed separately for each cohort's variant list", {
  long <- rbind(
    data.frame(cohort = "a", var_site = c("v1", "v2"),
               delta_wi_int = c(1, 2), log_odds_ratio = c(1, 2)),
    data.frame(cohort = "b", var_site = c("v1", "v2"),
               delta_wi_int = c(100, 200), log_odds_ratio = c(100, 200)))
  vs <- pgvis_variant_score(long)
  expect_equal(vs$variant_score[vs$cohort == "a"], vs$variant_score[vs$cohort == "b"])
})

test_that("delta AF is not in the default score, and requesting it warns", {
  long <- data.frame(cohort = "a", var_site = c("v1", "v2", "v3"),
                     delta_wi_int = 1:3, log_odds_ratio = 1:3, delta_af_alt = c(9, 1, 5))
  vs <- pgvis_variant_score(long)
  expect_false("p_delta_af_alt" %in% names(vs))
  expect_message(pgvis_variant_score(long, features = c("delta_wi_int", "delta_af_alt")),
                 "circular")
})

test_that("the sample score is the root mean square of the variant scores it carries", {
  long <- data.frame(sample = c("S1", "S1", "S2"), cohort = "a",
                     var_site = c("v1", "v2", "v1"),
                     variant_score = c(0.6, 0.8, 0.5))
  ps <- pgvis_sample_scores(long, by = c("sample", "cohort"))
  expect_equal(ps$pgvis_raw[ps$sample == "S1"], sqrt((0.36 + 0.64) / 2))
  expect_equal(ps$pgvis_raw[ps$sample == "S2"], 0.5)
  expect_equal(ps$n_variants, c(2L, 1L))
})

test_that("a variant listed twice for one sample is counted once", {
  long <- data.frame(sample = "S1", cohort = "a", var_site = c("v1", "v1", "v2"),
                     variant_score = c(0.6, 0.6, 0.8))
  ps <- pgvis_sample_scores(long, by = c("sample", "cohort"))
  expect_equal(ps$n_variants, 2L)
})

test_that("the sample score does not follow the number of variants carried", {
  set.seed(2)
  lg <- toy_long(enrich = FALSE)
  ps <- pgvis_sample_scores(pgvis_variant_score(lg))
  expect_lt(abs(stats::cor(ps$pgvis_raw, ps$n_variants, method = "spearman")), 0.5)
})

test_that("the reference scores near zero within ancestry and the test cohort sits above it", {
  set.seed(3)
  ps <- pgvis_sample_scores(pgvis_variant_score(toy_long()))
  s <- pgvis_score(ps, reference_group = "1000Genomes")
  ref <- s[s$cohort == "1000Genomes", ]
  tum <- s[s$cohort == "tumour", ]
  # reference is standardised against itself within superpopulation
  expect_equal(as.vector(tapply(ref$pgvis, ref$superpopulation, mean)), c(0, 0),
               tolerance = 1e-8)
  expect_equal(as.vector(tapply(ref$pgvis, ref$superpopulation, stats::sd)), c(1, 1),
               tolerance = 1e-8)
  expect_gt(median(tum$pgvis), 0.5)
  expect_equal(attr(s, "method"), "percentile")
  expect_true(all(c("mu", "sigma") %in% names(attr(s, "reference"))))
})

test_that("a frozen reference reproduces the scores", {
  set.seed(4)
  ps <- pgvis_sample_scores(pgvis_variant_score(toy_long()))
  s1 <- pgvis_score(ps)
  s2 <- pgvis_score(ps, reference = attr(s1, "reference"))
  expect_equal(s1$pgvis, s2$pgvis)
})

test_that("percentile scoring asks for pgvis_raw and names the way to get it", {
  d <- toy()
  expect_error(pgvis_score(d), "pgvis_sample_scores")
})

test_that("scoring fails loudly without a reference group", {
  set.seed(5)
  ps <- pgvis_sample_scores(pgvis_variant_score(toy_long()))
  ps$cohort <- "tumour"
  expect_error(pgvis_score(ps, reference_group = "1000Genomes"), "not present")
})

test_that("strata = NULL disables ancestry matching", {
  set.seed(6)
  ps <- pgvis_sample_scores(pgvis_variant_score(toy_long()))
  s <- pgvis_score(ps, strata = NULL)
  expect_true("pgvis" %in% names(s))
  expect_equal(nrow(s), nrow(ps))
})

# --- legacy PCA composite -----------------------------------------------------

test_that("legacy: the reference scores below the test cohort, and the pooled score is centred", {
  s <- pgvis_score(toy(), reference_group = "1000Genomes", method = "pca")
  ref <- s$pgvis[s$cohort == "1000Genomes"]
  tum <- s$pgvis[s$cohort == "tumour"]
  expect_lt(median(ref), median(tum))
  expect_gt(median(tum) - median(ref), 0.8)
  expect_lt(abs(mean(s$pgvis)), 0.1)
})

test_that("legacy: weights are PC1 loadings, non-negative and sum to one", {
  s <- pgvis_score(toy(), method = "pca")
  w <- attr(s, "weights")
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_true(all(w >= 0))
  expect_length(w, 3L)
})

test_that("legacy: frozen weights can be reused on a new cohort", {
  w <- attr(pgvis_score(toy(), method = "pca"), "weights")
  s2 <- pgvis_score(toy(), weights = w, method = "pca")
  expect_equal(attr(s2, "weights"), w)
})

test_that("carriers explode into one row per sample", {
  df <- data.frame(var_site = c("v1", "v2"),
                   carrier_samples = c("S1;S2,S3", "."),
                   stringsAsFactors = FALSE)
  out <- pgvis_explode_carriers(df)
  expect_equal(nrow(out), 3L)
  expect_setequal(out$sample, c("S1", "S2", "S3"))
})

test_that("superpopulation mapping leaves unknown values as NA", {
  m <- pgvis_map_superpopulation(c("white", "asian", "not reported", "other"))
  expect_equal(as.character(m), c("EUR", "EAS", NA, NA))
})

test_that("group comparison suppresses statistics below the minimum n", {
  d <- data.frame(pgvis = rnorm(20), group = rep(c("a", "b"), 10))
  r <- pgvis_compare_groups(d, min_n = 50)
  expect_false(r$valid)
  expect_true(is.na(r$cohens_d))
})

test_that("group comparison detects a real difference", {
  d <- data.frame(pgvis = c(rnorm(100, 0), rnorm(100, 2)),
                  group = rep(c("a", "b"), each = 100))
  r <- pgvis_compare_groups(d, min_n = 50, n_perm = 200)
  expect_true(r$valid)
  expect_gt(r$cohens_d, 1)
  expect_lt(r$perm_p, 0.05)
})
