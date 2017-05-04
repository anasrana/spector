#' Create regularly spaced regions along the whole genome.
#'
#' @param genome_version character. The genome version of your bam file.
#'        Currently supported are `"hg19"` and `"hg38"`.
#' @param region_size integer. Size of the regions computed along the genome.
#'        Should ideally be a power of 2, if it isn't it is coerced to one.
#' @param reg_overlap numeric. fraction of overlap for regions.
#'
#' @return A `tbl_df` with chromosome, start, and end columns.
#'
#' @export
#'
#' @importFrom dplyr filter mutate
#' @importFrom purrr map2
#' @importFrom tidyr unnest
#' @importFrom stringr str_c
#'
#' @family data import functions
#'
#' @examples
#' full_genome_regions("hg19", region_size = 2^20)
#'
full_genome_regions <- function(genome_version, region_size = NULL,
                            reg_overlap = 0) {

  if (is.null(region_size)) {
    message("region_size not specified, choosing default value 2^16")
    region_size <- 2^16
  } else if (region_size < 2^10) {
    stop("Region size below 2^10 is not supported,",
         " the LAS values are not meaningful.", call. = FALSE)
  } else {
    region_size <- 2^(floor(log2(min(region_size))))
  }

  if (reg_overlap > 1 | reg_overlap < 0) {
    stop("reg_overlap = ", reg_overlap,
         ": It should be a fraction of the region size, in the range [0, 1]",
         call. = FALSE)
  }

  if (!(genome_version %in% unique(genome_gap$genome))) {
    stop("Currently only `", str_c(unique(genome_gap$genome), collapse = ", "),
         "` are supported as genome_version", call. = FALSE )
  }

  region_df <-
  genome_size %>%
    filter(genome == genome_version) %>%
    mutate(gap = lapply(chrom, function(ch) {
            genome_gap %>%
              filter(chrom == ch & genome == genome_version)}),
           regions_chr = regionCompute(size_chrom = size,
                                       region_size = region_size,
                                       overlap_s = reg_overlap),
           regions_chr = map2(regions_chr, gap, filterGap)) %>%
    select(chrom, regions_chr) %>%
    unnest()

  return(region_df)
}


#' @importFrom tidyr separate_rows
#' @importFrom dplyr mutate group_by select summarise
#' @importFrom stringr str_c
#'
filterGap <- function(region, gap) {

  region %>%
    mutate(gap_st = str_c(gap$chromStart, collapse = ","),
           gap_ed = str_c(gap$chromEnd, collapse = ",")) %>%
    separate_rows(gap_st, gap_ed, sep = ",", convert = TRUE) %>%
    mutate(drop = (start >= gap_st & start <= gap_ed) |
                  (gap_st >= start & gap_st <= end) |
                  (gap_ed >= start & gap_st <= end)) %>%
    group_by(start, end) %>%
    summarise(drop = any(drop)) %>%
    filter(!drop) %>%
    group_by() %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    select(start, end)
}

#' @importFrom tibble data_frame
#'
regionCompute <- function(size_chrom = NULL, region_size = 2^14,
                          overlap_s = 0, gap_df = NULL) {
  lapply(size_chrom, function(i_size) {
    idx_n <- seq.int(from = 1, to = floor(i_size / region_size), by = 1)
    data_frame(
      start = (idx_n - 1) * region_size + 1 -
              (idx_n - 1) * ( overlap_s * region_size),
      end = idx_n * region_size -
            (idx_n - 1) * (overlap_s * region_size)
      )
  })
}
