context("QC")

s1_path <- spector_sample("sample1.bam")
s2_path <- spector_sample("sample2.bam")
id_path <- spector_sample("sample_id.txt")
basic_path <- spector_sample("basic.bed")

results_df <- readr::read_csv("result_basic-bed.csv", col_types = "cdidc")

test_that("spector complete run with giab", {
  results_test <- spector_qc(f_bam = spector_sample(""), region_size = 2^14)

  expect_that(nrow(results_test), equals(156))
  expect_that(ncol(results_test), equals(5))
})

test_that("spector verify results of complete run with custom bed file", {
  results_test <- spector_qc(f_bam = spector_sample(""), f_bed = basic_path)

  expect_equal(results_test, results_df)
})

test_that("spector verify results when using a parameter file", {
  results_test <- spector_qc(f_bam = id_path, f_bed = basic_path,
    file_type = "list")

  expect_equal(results_test[, 1:5], results_df)
})

test_that("chrIntersect check if intersect works in mixed cases", {
  bed_df <- read_bed(basic_path)

  mixed_bam <- nReadsBam(s2_path)
  basic_bam <- nReadsBam(s1_path)

  region_mix <- chrIntersect(bed_df, mixed_bam$chrom)
  region_smpl <- chrIntersect(bed_df, basic_bam$chrom)

  expect_that(unique(region_mix$chrom), equals(c("3", "4", "chr17")))
  expect_that(unique(region_smpl$chrom), equals(c("1", "2", "3", "4")))
  expect_that(nrow(region_mix), equals(12))
  expect_that(nrow(region_smpl), equals(17))
})
