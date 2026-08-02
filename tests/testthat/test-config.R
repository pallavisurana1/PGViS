test_that("config builds a full directory map and validates the genome", {
  cfg <- pgvis_config(root = file.path(tempdir(), "proj1"))
  expect_s3_class(cfg, "pgvis_config")
  expect_true(all(c("bam", "gvcf", "score", "logs") %in% names(cfg$dirs)))
  expect_error(pgvis_config(root = tempdir(), genome = "hg19"), "hg38")
})

test_that("thresholds merge with defaults rather than replacing them", {
  cfg <- pgvis_config(root = tempdir(), thresholds = list(ref_prob_min = 0.7))
  expect_equal(cfg$thresholds$ref_prob_min, 0.7)
  expect_equal(cfg$thresholds$alt_prob_max, 0.5)   # default retained
  expect_equal(cfg$thresholds$max_indel_len, 10L)
})

test_that("pgvis_init creates the tree", {
  root <- file.path(tempdir(), paste0("proj", as.integer(runif(1, 1, 1e6))))
  cfg <- pgvis_init(pgvis_config(root = root))
  expect_true(dir.exists(cfg$dirs$score))
})

test_that("unnamed regions are rejected", {
  expect_error(pgvis_config(root = tempdir(), regions = list("a.bed")), "named")
})

test_that("schemas have the expected widths", {
  expect_length(pgvis_schema("splice", with_dbsnp = FALSE), 8 + 12)
  expect_length(pgvis_schema("splice", with_dbsnp = TRUE), 8 + 12 + 17)
  expect_length(pgvis_schema("tfbs", with_dbsnp = TRUE), 8 + 12 + 17)
  expect_equal(pgvis_schema("splice")[1], "site_chr")
  expect_equal(pgvis_schema("splice", region_first = FALSE)[1], "var_chr")
})

test_that("all bundled shell stages exist on disk", {
  st <- pgvis_stages()
  expect_gt(nrow(st), 5)
  expect_true(all(file.exists(st$path)))
  expect_false(any(is.na(st$description)))
})
