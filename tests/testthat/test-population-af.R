make_row <- function(freqs, minor, ref, alt) {
  data.frame(dbsnp_freqs = freqs, dbsnp_minor = minor,
             dbsnp_major = minor, dbsnp_ref = ref, dbsnp_alt = alt,
             stringsAsFactors = FALSE)
}

test_that("TOPMED is preferred over lower-priority sources", {
  # index 1 = 1000G, 2 = dbGAP, 3 = TOPMED, 13 = GnoMAD
  f <- paste(c("0.11", "0.12", "0.13", rep("", 9), "0.19"), collapse = ",")
  m <- paste(c("A", "A", "A", rep("", 9), "A"), collapse = ",")
  out <- pgvis_resolve_population_af(make_row(f, m, "G", "A"))
  expect_equal(out$freq_source, "TOPMED")
  expect_equal(out$pop_af_alt, 0.13)
  expect_equal(out$pop_af_ref, 0.87)
})

test_that("a minor allele matching REF assigns the frequency to the reference", {
  f <- paste(c("", "", "0.30", rep("", 10)), collapse = ",")
  m <- paste(c("", "", "G", rep("", 10)), collapse = ",")
  out <- pgvis_resolve_population_af(make_row(f, m, "G", "A"))
  expect_equal(out$pop_af_ref, 0.30)
  expect_equal(out$pop_af_alt, 0.70)
})

test_that("sources fall through when the minor allele matches neither allele", {
  f <- paste(c("", "", "0.30", rep("", 9), "0.25"), collapse = ",")
  m <- paste(c("", "", "T", rep("", 9), "A"), collapse = ",")
  out <- pgvis_resolve_population_af(make_row(f, m, "G", "A"))
  expect_equal(out$freq_source, "GnoMAD")
  expect_equal(out$pop_af_alt, 0.25)
})

test_that("no usable source yields NA rather than a fabricated frequency", {
  out <- pgvis_resolve_population_af(make_row(",,,,,,,,,,,,", ",,,,,,,,,,,,", "G", "A"))
  expect_equal(out$freq_source, "none")
  expect_true(is.na(out$pop_af_alt))
})

test_that("priority can be reordered", {
  f <- paste(c("0.11", "", "0.13", rep("", 10)), collapse = ",")
  m <- paste(c("A", "", "A", rep("", 10)), collapse = ",")
  out <- pgvis_resolve_population_af(make_row(f, m, "G", "A"),
                                     priority = c("1000G", "TOPMED"))
  expect_equal(out$freq_source, "1000G")
  expect_equal(out$pop_af_alt, 0.11)
})
