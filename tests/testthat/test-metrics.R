test_that("log_odds_ratio is positive for disruptive variants and antisymmetric", {
  expect_gt(log_odds_ratio(0.95, 0.05), 0)
  expect_lt(log_odds_ratio(0.05, 0.95), 0)
  expect_equal(log_odds_ratio(0.7, 0.7), 0)
  expect_equal(log_odds_ratio(0.9, 0.2), -log_odds_ratio(0.2, 0.9))
})

test_that("log_odds_ratio clamps saturated probabilities instead of returning Inf", {
  expect_true(is.finite(log_odds_ratio(1, 0)))
  expect_true(is.finite(log_odds_ratio(0, 1)))
})

test_that("score_change is negative for loss of function and down-weights low probabilities", {
  expect_lt(score_change(0.9, 0.1), 0)
  expect_gt(score_change(0.1, 0.9), 0)
  # same absolute probability difference, but lower confidence overall
  expect_lt(abs(score_change(0.3, 0.1)), abs(score_change(0.9, 0.7)))
})

test_that("rank_inverse_normal is monotone, centred and NA-preserving", {
  x <- c(1, 5, 2, 100, 3)
  z <- rank_inverse_normal(x)
  expect_equal(order(z), order(x))
  expect_lt(abs(mean(z)), 1e-8)
  y <- rank_inverse_normal(c(1, NA, 3, 2))
  expect_true(is.na(y[2]))
  expect_equal(sum(!is.na(y)), 3L)
})

test_that("delta_wi is zero when cohort matches population and grows with divergence", {
  same <- data.frame(af_ref = 0.8, af_alt = 0.2,
                     pop_af_ref = 0.8, pop_af_alt = 0.2)
  expect_lt(pgvis_delta_wi(same, transform = FALSE, keep_raw = TRUE)$delta_wi, 1e-4)

  d <- data.frame(af_ref = c(0.8, 0.5), af_alt = c(0.2, 0.5),
                  pop_af_ref = c(0.8, 0.95), pop_af_alt = c(0.2, 0.05))
  raw <- pgvis_delta_wi(d, transform = FALSE, keep_raw = TRUE)$delta_wi
  expect_gt(raw[2], raw[1])
})

test_that("delta_wi derives the missing population reference frequency", {
  d <- data.frame(af_ref = 0.6, af_alt = 0.4, pop_af_alt = 0.1)
  out <- pgvis_delta_wi(d, transform = FALSE)
  expect_equal(out$pop_af_ref, 0.9)
  expect_equal(out$delta_af_alt, 0.3)
})

test_that("pgvis_filter_disruptive applies thresholds and records an audit trail", {
  df <- data.frame(
    p_ref = c(0.9, 0.9, 0.3, 0.9),
    p_alt = c(0.1, 0.8, 0.1, 0.1),
    n_carrier = c(100L, 100L, 100L, 1L),
    freq_source = c("TOPMED", "TOPMED", "TOPMED", "TOPMED")
  )
  out <- pgvis_filter_disruptive(df, ref_prob_min = 0.5, alt_prob_max = 0.5,
                                 min_carrier_frac = 0.05)
  expect_equal(nrow(out), 1L)          # only row 1 passes everything
  audit <- attr(out, "pgvis_filter")
  expect_equal(unname(audit["input"]), 4L)
  expect_equal(unname(audit["after_model_prob"]), 2L)
})
