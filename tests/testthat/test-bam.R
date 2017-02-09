context("Reading bam file")

s1_path <- spector_sample("sample1.bam")
s2_path <- spector_sample("sample2.bam")

test_that("", {
  s1_reads <- nReadsBam(s1_path)
  s2_reads <- nReadsBam(s2_path)

  expect_equal(nrow(s1_reads), 4)
  expect_equal(nrow(s2_reads), 3)
  expect_equal(s1_reads$n_read, c(111, 131, 152, 113))
  expect_equal(s2_reads$n_read, c(8, 157, 122))
})

test_that("chrCov reading coverage per chromosome", {
  res_cov <- read.csv("sample1_chr1_cov.csv", stringsAsFactor = F)

  cov <- spector:::chrCov(s1_path, "1", c(12762966, 47374980, 50387057, 241734214),
    c(12785179, 47387403, 50400445, 241747105), 1) %>%
  lapply(function(x) as.vector(x))

  cov_x <- strsplit(res_cov$x, ",") %>%
                    lapply(function(x) as.numeric(x))

  expect_that(cov, equals(cov_x))
})
