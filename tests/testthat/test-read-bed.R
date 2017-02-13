context("Reading bed file")

res_path <- spector_sample("bed_results.bed")
basic_path <- spector_sample("basic.bed")
header_path <- spector_sample("with_header.bed")

bed_df <- readr::read_csv(res_path, col_types = "cii")

test_that("read_bed split size, correct region split", {
  bed_test <- read_bed(basic_path)

  expect_that(nrow(bed_test), equals(21))
  expect_that(bed_test$chrom, equals(bed_df$chrom))
  expect_true(is.numeric(bed_test$start))
  expect_true(is.numeric(bed_test$end))
  expect_that(bed_test, equals(bed_df))
})

test_that("read_bed including header", {
  bed_test <- read_bed(header_path, header = TRUE)

  expect_that(nrow(bed_test), equals(21))
  expect_that(bed_test, equals(bed_df))
})

test_that("read_bed ucsc coordinate shift", {
  bed_test <- read_bed(header_path, header = TRUE, ucsc_coord = TRUE)

  expect_that(nrow(bed_test), equals(21))
  expect_that(bed_test$start[1], equals(12764911))
  expect_that(bed_test$start[17], equals(133006534))
  expect_that(bed_test$end[18], equals(14402))
  expect_that(bed_test$end[13], equals(179884855))
  expect_is(bed_test$start, "integer")
  expect_is(bed_test$end, "integer")
  expect_is(bed_test$chrom, "character")
})

test_that("read_bed custom regions_size", {
  bed_test <- read_bed(header_path, header = TRUE, ucsc_coord = TRUE,
    bed_region_size = 2^10)
  expect_that(nrow(bed_test), equals(236))
  expect_that(ncol(bed_test), equals(3))
})

test_that("checkRegionSize test if correct region supplied", {
  region_size <- checkRegionSize(2000,
                  dplyr::data_frame(reg_length = c(12, 20, 1999, 2000, 5000)))

  expect_that(region_size, equals(1024))
})

test_that("full genome regions computed accuratley", {
  region_df <- full_genome_regions(genome_version = "hg19", region_size = 2^20)

  expect_equal(nrow(region_df), 2628)
  expect_equal(ncol(region_df), 3)
  expect_is(region_df$chrom, "character")
  expect_is(region_df$start, "integer")
  expect_is(region_df$end, "integer")
  expect_error(full_genome_regions(genome_version = "hg18"))
  expect_error(full_genome_regions(genome_version = "hg19",
                                   region_size = "2^9"))
})
