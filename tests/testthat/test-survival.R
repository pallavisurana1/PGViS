test_that("stratification produces the requested number of groups", {
  x <- rnorm(90)
  expect_length(levels(pgvis_stratify(x, "median")), 2L)
  expect_length(levels(pgvis_stratify(x, "tertile")), 3L)
  expect_length(levels(pgvis_stratify(x, "quartile")), 4L)
  expect_equal(as.character(pgvis_stratify(c(1, 3), "cutpoint", cutpoint = 2)),
               c("low", "high"))
})

test_that("stratification refuses to split heavily tied data", {
  expect_error(pgvis_stratify(rep(1, 50), "tertile"), "tied")
})

test_that("survival runs and detects a planted effect", {
  skip_if_not_installed("survival")
  set.seed(3)
  n <- 300
  score <- rnorm(n)
  # higher score -> shorter time
  time <- rexp(n, rate = exp(0.7 * score) / 500)
  scored <- data.frame(sample = paste0("S", 1:n), pgvis = score)
  clin <- data.frame(sample = paste0("S", 1:n), time = time,
                     event = rbinom(n, 1, 0.8))
  res <- pgvis_survival(scored, clin)
  expect_true(all(c("fit", "cox_continuous", "summary") %in% names(res)))
  expect_lt(res$summary$cox_p, 0.01)
  expect_gt(res$summary$hr, 1)
})

test_that("survival rejects non-binary event coding", {
  skip_if_not_installed("survival")
  scored <- data.frame(sample = "S1", pgvis = 1)
  clin <- data.frame(sample = "S1", time = 10, event = 2)
  expect_error(pgvis_survival(scored, clin), "0/1")
})
