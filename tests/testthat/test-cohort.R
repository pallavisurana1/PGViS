test_that("sample sheets round-trip and reject breakage", {
  x <- data.frame(cohort = c("luad", "lusc"), tissue = "blood",
                  sample = c("S1", "S2"))
  p <- file.path(tempdir(), "ss.tsv")
  pgvis_sample_sheet(x, path = p)
  back <- pgvis_sample_sheet(p)
  expect_equal(nrow(back), 2L)
  expect_equal(back$sample, c("S1", "S2"))

  expect_error(pgvis_sample_sheet(data.frame(cohort = "a", sample = "b")),
               "Missing")
  expect_error(pgvis_sample_sheet(
    data.frame(cohort = "a\tb", tissue = "t", sample = "s")), "Whitespace")
  expect_error(pgvis_sample_sheet(
    data.frame(cohort = "a", tissue = "t", sample = c("s", "s"))), "duplicated")
})

test_that("composite labels split into components without cross-contamination", {
  p <- pgvis_parse_label(c("tcga_luad_smoke", "cptac_lusc", "1000Genomes"))
  expect_equal(p$study,   c("tcga", "cptac", NA))
  expect_equal(p$cancer,  c("luad", "lusc", NA))
  expect_equal(p$stratum, c("smoke", NA, NA))
})

test_that("replication filter keeps only sites seen in two studies", {
  df <- data.frame(
    study     = c("tcga", "cptac", "tcga", "tcga"),
    cancer    = "luad",
    func_site = c("f1", "f1", "f2", "f2"),
    var_site  = c("v1", "v1", "v2", "v2"),
    stringsAsFactors = FALSE
  )
  out <- pgvis_filter_replicated(df)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$func_site == "f1"))
})

test_that("replication is required within stratum, not across it", {
  df <- data.frame(study = c("tcga", "cptac"), cancer = c("luad", "lusc"),
                   func_site = "f1", var_site = "v1", stringsAsFactors = FALSE)
  expect_equal(nrow(pgvis_filter_replicated(df, within = "cancer")), 0L)
  expect_equal(nrow(pgvis_filter_replicated(df, within = character())), 2L)
})
