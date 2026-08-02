test_that("QC flags degenerate model predictions", {
  df <- data.frame(p_ref = rep(0.5, 100), p_alt = rep(0.5, 100))
  qc <- pgvis_qc(df, verbose = FALSE)
  expect_true(qc$degenerate_predictions)
  expect_true(any(grepl("distinct values", qc$warnings)))
})

test_that("QC flags mostly-missing population frequencies", {
  df <- data.frame(freq_source = c(rep("none", 80), rep("TOPMED", 20)))
  qc <- pgvis_qc(df, verbose = FALSE)
  expect_equal(qc$no_population_af, 80L)
  expect_true(any(grepl("population frequency", qc$warnings)))
})

test_that("QC on the clean example produces no warnings", {
  v <- pgvis_example("variants"); p <- pgvis_example("predictions")
  m <- pgvis_variant_metrics(v, p)
  qc <- pgvis_qc(m, verbose = FALSE)
  expect_length(qc$warnings, 0L)
})

test_that("the manifest records thresholds", {
  cfg <- pgvis_init(pgvis_config(root = file.path(tempdir(), "man-test"),
                                 thresholds = list(ref_prob_min = 0.7)))
  m <- pgvis_manifest(cfg, "test", inputs = list(n = 10))
  expect_equal(m$thresholds$ref_prob_min, 0.7)
  expect_true(length(list.files(cfg$dirs$logs, pattern = "^manifest_")) >= 1)
})
