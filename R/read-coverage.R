#' Read the coverage for regions on one chromosome.
#'
#' @param f_name string. Full path to bam file.
#' @param chr string. Chromosome name, make sure it follows to naming
#'        standards in your bam file.
#' @param start integer vector. Starting coordinates for regions to read.
#' @param end integer vector. Starting coordinates for regions to read.
#' @param n_read integer. total number of reads in the bam file.
#'
#' @return GRanges object with reads from a given chromosome.
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignments coverage
#' @importFrom magrittr extract2
#'
#' @family data import functions
#'
#' @export
#'
read_cov <- function(f_name, chr, start, end, n_read = NULL) {
  # Create param from provided information
  param <- ScanBamParam(which = rngObj(chr, start, end))

  raw_reads <- readGAlignments(f_name, param = param,
      use.names = F) %>%
    coverage() %>%
    extract2(chr)

  raw_reads <- raw_reads / n_read
  return(raw_reads)
}

#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#'
rngObj <- function(chr, start, end) {
  GRanges(Rle(chr),
    ranges = IRanges(start = as.integer(start), end = as.integer(end)))
}
