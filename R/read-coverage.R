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
read_cov <- function(f.name, chr, start, end, n.read) {
  # Create param from provided information
  flag <- Rsamtools::scanBamFlag(isNotPassingQualityControls=FALSE,
    isDuplicate=FALSE, isProperPair = TRUE, hasUnmappedMate = FALSE)

  param <- Rsamtools::ScanBamParam(which = rng_obj(chr, start, end),
    flag = flag, mapqFilter = 30)

  raw.reads <- GenomicAlignments::readGAlignments(f.name, param = param,
      use.names = F) %>%
    coverage() %>%
    magrittr::extract2(chr) %>%
    as.numeric() %>%
    magrittr::extract(seq(start, end, by = 1))

  read.norm <- raw.reads / n.read

  if (log2(length(read.norm)) == round(log2(length(read.norm)))) {
        signal <- read.norm
  } else {
      if (length(read.norm) %% 2 == 0) {
          ll <- length(read.norm) / 2 - 2^floor(log2(length(read.norm))) / 2 + 1
          ul <- length(read.norm) / 2 + 2^floor(log2(length(read.norm))) / 2
      } else {
          ll <- round(length(read.norm) / 2) -
                2^floor(log2(length(read.norm))) / 2 + 1
          ul <- round(length(read.norm) / 2) + 2^floor(log2(length(read.norm))) / 2
      }
      signal <- read.norm[ll:ul]
  }

  signal <- .sig_norm(signal)

  return(signal)
}

.sig_norm <- function(signal) {
  if (max(signal, na.rm = TRUE) != 0) {
    signal <- (signal - min(signal, na.rm = TRUE)) /
              (max(signal, na.rm = TRUE) - min(signal, na.rm = TRUE))
  }

  return(signal)
}

rng_obj <- function(chr, start, end) {
  GenomicRanges::GRanges(Rle(chr),
    ranges = IRanges(start = as.integer(start), end = as.integer(end)))
}
