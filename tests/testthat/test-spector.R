context("QC")

results_df <- readr::read_csv("result_basic-bed.csv", col_types = "cdidc")

test_that("spector complete run with giab", {
  results_test <- spector_qc(f_bam = ".", region_size = 2^14)

  expect_that(nrow(results_test), equals(156))
  expect_that(ncol(results_test), equals(5))
})

test_that("spector verify results of complete run with custom bed file", {
  results_test <- spector_qc(f_bam = ".", f_bed = "basic.bed")

  expect_equal(results_test, results_df)
})

test_that("chrIntersect check if intersect works in mixed cases", {
  bed_df <- read_bed("basic.bed")

  mixed_bam <- nReadsBam("sample2.bam")
  basic_bam <- nReadsBam("sample1.bam")

  region_mix <- chrIntersect(bed_df, mixed_bam$chrom)
  region_smpl <- chrIntersect(bed_df, basic_bam$chrom)

  expect_that(unique(region_mix$chrom), equals(c("3", "4", "chr17")))
  expect_that(unique(region_smpl$chrom), equals(c("1", "2", "3", "4")))
  expect_that(nrow(region_mix), equals(12))
  expect_that(nrow(region_smpl), equals(17))
})
