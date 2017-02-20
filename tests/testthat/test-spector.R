context("Main QC function")

s1_path <- spector_sample("sample1.bam")
s2_path <- spector_sample("sample2.bam")
id_path <- spector_sample("sample_id.txt")
id_path_2 <- spector_sample("sample_id_2.txt")
basic_path <- spector_sample("basic.bed")

results_df <- readr::read_csv("result_basic-bed_rms.csv", col_types = "ciidcc")
results_all <- readr::read_csv("result_basic-bed.csv", col_types = "ciidddcc")

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

  expect_equal(results_test[, 1:6], results_df)
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


test_that("test passing bed file as tbl_df", {
  basic_df <- read_bed(basic_path)
  results_test <- spector_qc(f_bam = spector_sample(""), f_bed = basic_df)

  expect_equal(results_test, results_df)
})

test_that("spector complete run with giab", {
  results_test <- spector_qc(f_bam = spector_sample(""), region_size = 2^16)

  expect_that(nrow(results_test), equals(109))
  expect_that(ncol(results_test), equals(6))
})

test_that("spector verify results of complete run with custom bed file", {
  results_test <- spector_qc(f_bam = spector_sample(""), f_bed = basic_path)

  expect_equal(results_test, results_df)
})

test_that("test spector on single file", {
  results_test <- spector_qc(f_bam = s1_path, f_bed = basic_path)

  expect_equal(nrow(results_test), 17)
  expect_equal(ncol(results_test), 6)
  expect_equal(results_test$metric[16], 7555236.720667904)
})

test_that("test basic run of spectorFile", {

  results_s2 <-
  spector:::spectorFile(f_bam = s2_path, id_bam = "file_id", s_prep = "prep_id",
                        out_F = NULL, save_out = FALSE, chr_cores = 1,
                        regions = "giab", region_size = 2^16, f_bed = NULL,
                        header = FALSE)

  results_s2_bed <-
  spector:::spectorFile(f_bam = s2_path, id_bam = "file_id", s_prep = "prep_id",
                        out_F = NULL, save_out = FALSE, chr_cores = 1,
                        regions = "giab", region_size = 2^14,
                        f_bed = basic_path, header = FALSE)

  expect_equal(nrow(results_s2), 35)
  expect_equal(ncol(results_s2), 7)
  expect_equal(unique(c(results_s2$id_bam, results_s2_bed$id_bam)), "file_id")
  expect_equal(unique(c(results_s2$prep, results_s2_bed$prep)), "prep_id")

})

test_that("spector runs with function returning all variables", {

  results_test <- spector_qc(f_bam = spector_sample(""), f_bed = basic_path,
                             smr_var = "all")

  expect_equal(results_test, results_all)

})

test_that("spector fails when region sizes small", {
  expect_error(spector_qc(s1_path, region_size = 2^9))
})

test_that("spector picks up genome version", {
expect_message(
   tmp <-
    spector:::spectorFile(f_bam = s2_path, id_bam = "file_id", s_prep = "prep_id",
                        out_F = NULL, save_out = FALSE, chr_cores = 1,
                        regions = "giab", region_size = 2^16, f_bed = NULL,
                        header = FALSE, genome_v = "hg38"),
  "Genome version selected: hg38")

expect_message(
  spector:::spectorFile(f_bam = s2_path, id_bam = "file_id", s_prep = "prep_id",
                      out_F = NULL, save_out = FALSE, chr_cores = 1,
                      regions = "giab", region_size = 2^16, f_bed = NULL,
                      header = FALSE, genome_v = "GRCh37"),
  "Genome version selected: hg19")

expect_message(
   spector_qc(f_bam = s2_path, f_bed = basic_path),
  "Genome version selected: hg19")

expect_message(
   spector_qc(f_bam = s2_path, f_bed = basic_path, genome = "GRCh38"),
  "Genome version selected: hg38")
})

test_that("spector fails with wrong genome", {

expect_error(
  spector:::spectorFile(f_bam = s2_path, id_bam = "file_id", s_prep = "prep_id",
                      out_F = NULL, save_out = FALSE, chr_cores = 1,
                      regions = "giab", region_size = 2^16, f_bed = NULL,
                      header = FALSE, genome_v = "hg18")
)

expect_error(spector_qc(f_bam = s2_path, f_bed = basic_path, genome = "hg18"))

})
