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

test_that("the reference cohort scores below the test cohort, and the pooled score is centred", {
  s <- pgvis_score(toy(), reference_group = "1000Genomes")
  ref <- s$pgvis[s$cohort == "1000Genomes"]
  tum <- s$pgvis[s$cohort == "tumour"]
  # The composite z-scores are re-scaled across the pooled cohort, so the
  # reference arm sits below zero whenever the test arm is shifted up; what
  # the score guarantees is separation and a centred pool, not a zero-centred
  # reference.
  expect_lt(median(ref), median(tum))
  expect_gt(median(tum) - median(ref), 0.8)
  expect_lt(abs(mean(s$pgvis)), 0.1)
})

test_that("weights are PC1 loadings, non-negative and sum to one", {
  s <- pgvis_score(toy())
  w <- attr(s, "weights")
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_true(all(w >= 0))
  expect_length(w, 3L)
})

test_that("frozen weights can be reused on a new cohort", {
  w <- attr(pgvis_score(toy()), "weights")
  s2 <- pgvis_score(toy(), weights = w)
  expect_equal(attr(s2, "weights"), w)
})

test_that("scoring fails loudly without a reference group", {
  d <- toy(); d$cohort <- "tumour"
  expect_error(pgvis_score(d, reference_group = "1000Genomes"), "not present")
})

test_that("strata = NULL disables ancestry matching", {
  s <- pgvis_score(toy(), strata = NULL)
  expect_true("pgvis" %in% names(s))
  expect_equal(nrow(s), 240L)
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
