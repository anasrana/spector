#' Read bam file with regions
#'
#' @param file Path to the bed file
#' @param header \code{TRUE} if the bed file contains a header
#' @param ucsc.coord \code{TRUE} if coordinates start with 0 like ucsc
#'
#' @return A \code{data.frame} used in the following calculations
#'
spector_bed <- function(file, header = FALSE, ucsc.coord = TRUE) {
# column names and first line
# ----------------------------------------------------------------

  c_name <- c("chrom", "start", "end", "name", "score", "shade", "strand",
   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
  tmp <- readr::read_delim(file, delim = "\t", n_max = 1, col_names = FALSE)

# Reading bed files
# ----------------------------------------------------------------
  region <- readr::read_delim(file, delim = "\t",  col_names = c_name[1:ncol(tmp)],
                              progress = FALSE,
                              col_type = list(
                                chrom = readr::col_character(),
                                start = readr::col_integer(),
                                end = readr::col_integer()
                                ))
# removing unneccsary columns
# ----------------------------------------------------------------
  if (ncol(region > 3)) {
    region <- region %>%
      dplyr::select(chrom, start, end)
  }

# Shift the data if ucsc convention
# ----------------------------------------------------------------

  if (is.integer(region$chrom)) {
    if (ucsc.coord) {
      region <- region %>%
        dplyr::mutate(start = start + 1,
          id = paste("chr", chrom, ":", start, "-", end, sep = ""),
          chr = paste("chr", chrom, sep = ""))
      } else {
        region <- region %>%
          dplyr::mutate(id = paste("chr", chrom, ":", start, "-", end, sep = ""),
            chr = paste("chr", chrom, sep = ""))
        }
    } else {
      if (ucsc.coord) {
        region <- region %>%
          dplyr::mutate(start = start + 1,
            id = paste("chr", chrom, ":", start, "-", end, sep = ""), chr = chrom)
        } else {
          region <- region %>%
            dplyr::mutate(id = paste("chr", chrom, ":", start, "-", end, sep = ""), chr = chrom)
        }
  }

  return(region)
}
