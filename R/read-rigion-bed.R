#' Read bam file with regions
#'
#' @param file Path to the bed file
#' @param header \code{TRUE} if the bed file contains a header
#' @param ucsc.coord \code{TRUE} if coordinates start with 0 like ucsc
#'
#' @return A \code{data.frame} used in the following calculations
#' @export
spector_bed <- function(file, header = FALSE, ucsc.coord = TRUE) {
# Reading data files ----------------------------------------------------------------
  region <- readr::read_delim(file, delim = "\t", col_names = header)
  colnames(region) <- c("chrom", "start", 'end')

# Shift the data if ucsc convention -------------------------------------------------
  if (is.integer(region$chrom)) {
    if (ucsc.coord) {
      region <- region %>%
        mutate(start = start + 1,
          id = paste("chr", chrom, ":", start, "-", end, sep = ""),
          chr = paste("chr", chrom, sep = ""))
      } else {
        region <- region %>%
          mutate(id = paste("chr", chrom, ":", start, "-", end, sep = ""),
            chr = paste("chr", chrom, sep = ""))
        }
    } else {
      if (ucsc.coord) {
        region <- region %>%
          mutate(start = start + 1,
            id = paste("chr", chrom, ":", start, "-", end, sep = ""), chr = chrom)
        } else {
          region <- region %>%
            mutate(id = paste("chr", chrom, ":", start, "-", end, sep = ""), chr = chrom)
        }
  }

  return(region)
}
