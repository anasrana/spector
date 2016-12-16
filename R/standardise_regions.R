#' Read bed bed_file with bed_regions to run spector
#'
#' @param bed_file Path to the bed file
#' @param header \code{TRUE} if the bed file contains a header
#' @param ucsc_coord \code{TRUE} if coordinates start with 0 like ucsc
#'
#' @return A \code{data.frame} used in the following calculations
#'
#' @importFrom dplyr as_data_frame select mutate
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#'
#' @export
#'
read_bed <- function(bed_file,
  header = FALSE,
  ucsc_coord = TRUE,
  bed_region_size = NULL) {

# column names and first line
# ----------------------------------------------------------------

  c_name <- c("chrom", "start", "end", "name", "score", "shade", "strand",
   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
   "blockStarts")

  if (header) {
    tmp.row1 <- read.table(bed_file , nrows = 1, skip = 1)
  } else {
    tmp.row1 <- read.table(bed_file , nrows = 1)
  }

# Reading bed files
# ----------------------------------------------------------------

  bed_region <- read.table(
      bed_file,
      as.is = TRUE,
      stringsAsFactors = F,
      sep = "\t") %>%
    as_data_frame() %>%
    setNames(c_name[1:ncol(tmp.row1)]) %>%
    select(chrom, start, end) %>%
    mutate(reg_length = end - start)

  if (is.null(bed_region_size)) {
    message("No region size specified:
      Using largest power of 2 that fits into min(region) in the bed file")
    bed_region_size <- 2^(floor(log2(min(bed_region$reg_length))))
  } else {
    bed_region_size <- 2^(floor(log2(min(bed_region_size))))
    message(paste0("Regions standardised to length ", bed_region_size))
  }

  if (max(bed_region$reg_length) < bed_region_size) {
    stop(paste0("min region size invalid
      Choose a number smaller than: ", max(bed_region$reg_length)))
  }

  bed_region <-
    bed_region %>%
    bed_region_split(bed_region_size)

# Shift the data if ucsc convention
# ----------------------------------------------------------------

  if (is.integer(bed_region$chrom)) {
    if (ucsc_coord) {
      bed_region <- bed_region %>%
        dplyr::mutate(start = start + 1,
          id = paste("chr", chrom, ":", start, "-", end, sep = ""),
          chr = paste("chr", chrom, sep = ""))
      } else {
        bed_region <- bed_region %>%
          dplyr::mutate(
            id = paste("chr", chrom, ":", start, "-", end, sep = ""),
            chr = paste("chr", chrom, sep = ""))
        }
    } else {
      if (ucsc_coord) {
        bed_region <- bed_region %>%
          dplyr::mutate(start = start + 1,
            id = paste("chr", chrom, ":", start, "-", end, sep = ""),
              chr = chrom)
        } else {
          bed_region <- bed_region %>%
            dplyr::mutate(
              id = paste("chr", chrom, ":", start, "-", end, sep = ""),
              chr = chrom)
        }
  }

  return(bed_region)
}


# ==========================================================================
# Determine number of bed_regions to be split into and uncovered area
# ==========================================================================

#' mutate function bed_region_split
#'
#' @param bed_region
#' @param bed_region_size
#'
#' @return
#'
#' @examples
#'
#' @importFrom dplyr as_data_frame select mutate filter
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#' @importFrom tidyr separate_rows separate
#'
bed_region_split <- function(bed_region, bed_region_size) {
  bed_region %>%
    filter(reg_length > bed_region_size) %>%
    mutate(
      n_reg = reg_length %/% bed_region_size,
      uncov = reg_length %% bed_region_size,
      new.reg = bed_expand_bed_region(bed_region_size, start, end, n_reg, uncov),
      orig_reg = str_c(chrom, ":", start, "-", end)
      ) %>%
    select(-start, -end, -reg_length) %>%
    separate_rows(new.reg, sep = ",") %>%
    separate(new.reg, into = c("start", "end"), convert = TRUE) %>%
    select(chrom, start, end)
}

# ==========================================================================
#
# ==========================================================================

#' Title
#'
#' @param bed_region_size
#' @param start
#' @param end
#' @param n_reg
#' @param uncov
#'
#' @return
#'
#' @examples
#'
#' @importFrom dplyr as_data_frame select mutate
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#'
bed_expand_bed_region <- function(bed_region_size, start, end, n_reg, uncov) {
sapply(1:length(n_reg), function(i_v) {
  gap.bp <- round(uncov[i_v] / (n_reg[i_v] + 1))
  tmp.st <- c(rep(NA, n_reg[i_v]))
  tmp.ed <- c(start[i_v], rep(NA, n_reg[i_v]))
  res.v <- rep(NA, n_reg[i_v])

  for (i_n in 1:n_reg[i_v]) {
    tmp.st[i_n] <- tmp.ed[i_n] + gap.bp
    tmp.ed[i_n + 1] <- tmp.st[i_n] + bed_region_size
    res.v[i_n] <- str_c(tmp.st[i_n], "-", tmp.ed[i_n + 1])
  }

  res.v %>% str_c(collapse = ",")
})

}
