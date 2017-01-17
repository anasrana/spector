#' Read the coverage of a region
#'
#' @param chr
#' @param start
#' @param end
#' @param n_read
#' @param f_name
#'
#' @return Vector of the coverage in a given region
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignments coverage
#' @importFrom magrittr extract2
#'
#' @export
#'
read_cov <- function(f_name, chr, start, end, n_read) {
  # Create param from provided information
  param <- ScanBamParam(which = rngObj(chr, start, end))

  raw_reads <- readGAlignments(f_name, param = param,
      use.names = F) %>%
    coverage() %>%
    extract2(chr)

    return(raw_reads)
}

#' generate genomic ranges object
#'
#' @param chr
#' @param start
#' @param end
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#'
rngObj <- function(chr, start, end) {
  GRanges(Rle(chr),
    ranges = IRanges(start = as.integer(start), end = as.integer(end)))
}
