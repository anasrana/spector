#' @importFrom stringr str_detect
#'
checkGenome <- function(genome) {
#ToDo check for implementation to allow for different case and mixed case
  if (str_detect(genome, "38|19|37")) {
    if (str_detect(genome, "38")) {
      genome <- "hg38"
    } else if (str_detect(genome, "37|19")) {
      genome <- "hg19"
    }
  }

  if (genome %in% c("hg19", "hg38")) {
    message("Genome version selected: ", genome)
  } else {
    stop("This package only supports `GRCh37` and `GRCh38` genomes.\n",
       "  For any other genome version please provide a relevant bed file ",
       "using the `f_bed =` option")
  }

}

#' @importFrom dplyr filter select summarise
#' @importFrom stringr str_c
#'
lasFilter <- function(r_size, coff = 0.8) {
  if (r_size < min(las_limits$region_size) |
      r_size > max(las_limits$region_size)) {
    stop("Supplied region size not supported for")
  }

  min_c <- las_limits %>%
    filter(region_size == r_size) %>%
    summarise(min_c = min(c_off)) %>%
    with(min_c)

  if (coff < min_c) {
    stop(str_c("Your cutoff need to be > ", min_c, " for this region_size."))
  }

  las_limits  %>%
    filter(region_size == r_size & c_off == coff ) %>%
    select(lower_lim, upper_lim)
}

regionSize <- function(reg_DF) {
  reg_DF$end[1] - reg_DF$start[1] + 1
}

lasLimBin <- function(las, las_lim) {

  las < las_lim$lower_lim | las > las_lim$upper_lim
}
