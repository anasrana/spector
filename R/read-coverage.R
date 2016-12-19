#' Read the coverage of a region
#'
#' @param chr
#' @param start
#' @param end
#' @param n.read
#' @param f.name
#'
#' @return Vector of the coverage in a given region
#' @export
#'
read_cov <- function(f.name, chr, start, end, n_read) {
  # Create param from provided information
  param <- Rsamtools::ScanBamParam(which = rng_obj(chr, start, end))

  raw_reads <- GenomicAlignments::readGAlignments(f.name, param = param,
      use.names = F) %>%
    IRanges::coverage() %>%
    magrittr::extract2(chr)

    return(raw_reads)

  # read.norm <- raw.reads / n.read

  # if (log2(length(read.norm)) == round(log2(length(read.norm)))) {
  #       signal <- read.norm
  # } else {
  #     if (length(read.norm) %% 2 == 0) {
  #         ll <- length(read.norm) / 2 - 2^floor(log2(length(read.norm))) / 2 + 1
  #         ul <- length(read.norm) / 2 + 2^floor(log2(length(read.norm))) / 2
  #     } else {
  #         ll <- round(length(read.norm) / 2) -
  #               2^floor(log2(length(read.norm))) / 2 + 1
  #         ul <- round(length(read.norm) / 2) +
  #             2^floor(log2(length(read.norm))) / 2
  #     }
  #     signal <- read.norm[ll:ul]
  # }

  # signal <- .sig_norm(signal)
  # signal <- .sig_norm(read.norm)

  # return(signal)
}

rng_obj <- function(chr, start, end) {
  GenomicRanges::GRanges(S4Vectors::Rle(chr),
    ranges = IRanges::IRanges(start = as.integer(start), end = as.integer(end)))
}
