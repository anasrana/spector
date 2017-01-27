context("Main QC function")

s1_path <- spector_sample("sample1.bam")
s2_path <- spector_sample("sample2.bam")
id_path <- spector_sample("sample_id.txt")
id_path_2 <- spector_sample("sample_id_2.txt")
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
  import_par <- read_par_file(id_path)
  spector:::unpackList(import_par)

  expect_true(exists("fs_bam") & !is.null(fs_bam))
  expect_true(exists("id_bam"))
  expect_true(exists("sample_type"))
  expect_equal(nrow(import_par), 2)

  fs_bam <- c(s1_path, s2_path)

  results_test <-
    spector:::spectorList(fs_bam = fs_bam, id_v = id_bam, s_v = sample_type,
                              out_F = NULL, file_cores = 1, chr_cores = 1,
                              save_out = FALSE, f_bed = basic_path)

  expect_equal(results_test[, 1:5], results_df)
  expect_equal(unique(results_test$prep), sample_type)
})

test_that("test parameter file with and without header", {
  import_par <- read_par_file(id_path)
  import_par2 <- read_par_file(id_path_2)

  expect_equal(ncol(import_par), ncol(import_par2))
  expect_equal(nrow(import_par), nrow(import_par2))
  expect_equal(import_par2[, -1], import_par[, -1])
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


test_that("test bassing bed file as tbl_df", {
  basic_df <- read_bed(basic_path)
  results_test <- spector_qc(f_bam = spector_sample(""), f_bed = basic_df)

  expect_equal(results_test, results_df)
})
