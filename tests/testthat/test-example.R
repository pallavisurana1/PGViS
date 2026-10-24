test_that("bundled example data is present and self-consistent", {
  v <- pgvis_example("variants")
  p <- pgvis_example("predictions")
  expect_gt(nrow(v), 100)
  expect_equal(nrow(v), nrow(p))
  expect_setequal(unique(v$feature), unique(p$feature))
  expect_true(all(c("func_site", "var_site", "af_alt", "pop_af_alt",
                    "carrier_samples") %in% names(v)))
  expect_true(all(c("p_ref", "p_alt") %in% names(p)))
  expect_true(all(p$p_ref > 0 & p$p_ref < 1))
})

test_that("clinical and pedigree fixtures load", {
  cl <- pgvis_example("clinical")
  expect_true(all(c("sample", "race", "time", "event") %in% names(cl)))
  expect_true(all(cl$event %in% c(0, 1)))

  ped <- pgvis_load_pedigree(
    system.file("extdata", "toy_pedigree.txt", package = "PGViS"))
  expect_true(all(c("sample", "superpopulation") %in% names(ped)))
  expect_gt(nrow(ped), 10)
})

test_that("the example workflow runs end to end and produces finite scores", {
  res <- pgvis_example_run(verbose = FALSE)
  expect_true(all(c("metrics", "per_sample", "scores") %in% names(res)))
  expect_gt(nrow(res$scores), 0)
  expect_true("pgvis" %in% names(res$scores))
  expect_true(any(is.finite(res$scores$pgvis)))
  expect_true("pgvis_raw" %in% names(res$scores))
  expect_true(all(res$scores$pgvis_raw[!is.na(res$scores$pgvis_raw)] > 0 &
                  res$scores$pgvis_raw[!is.na(res$scores$pgvis_raw)] <= 1))
  ref <- attr(res$scores, "reference")
  expect_true(all(c("mu", "sigma") %in% names(ref)))
})

test_that("metrics computed on the example are in sensible ranges", {
  v <- pgvis_example("variants"); p <- pgvis_example("predictions")
  m <- pgvis_variant_metrics(v, p)
  expect_true(all(is.finite(m$log_odds_ratio)))
  expect_true(all(m$score_change <= 1 & m$score_change >= -1))
  expect_true(all(is.finite(stats::na.omit(m$delta_wi_int))))
})
