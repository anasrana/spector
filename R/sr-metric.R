#' @details
#' The \code{spector_metric} function is the main function that calculates a
#' metric value for each region in the selected bed file.
#'
#' The \code{spector_metric} function performs the actual calculations, but does
#' not include any checks or verifications so should ideally not be called
#' directly. It should be accessed via \code{\link{spector_qc}}.
#'
#' @param f_bam path to bam file, relative to current working directory
#'        (see \code{\link[base]{setwd}()} for more details).
#' @param region_size integer. Choose size of regions to calculate metric value.
#'        The default \code{ = NULL} means \code{region_size = } maximum power
#'        of 2 that fits in the smallest region.
#' @param f_bed file path for bed file to override default giab.
#' @param bed_header logical. \code{TRUE} if bed file has a header, the default
#'        value is \code{FALSE}.
#' @param region_giab logical. Indicates whether or not giab regions are used,
#'        defaults to \code{TRUE}.
#' @param chr_cores integer. Optional number indicating if the QC should be
#'        computed in parallel across chromosomes. Default value is \code{1}.
#'
#' @return Object (\code{tbl_df}) of metric for all specified regions.
#'
#' @importFrom dplyr mutate group_by summarise bind_rows filter
#' @importFrom stringr str_c
#' @importFrom magrittr %>%
#' @importFrom tidyr separate separate_rows
#' @importFrom parallel mclapply
#'
#' @export
#'
spector_metric <- function(region_df, f_bam = NULL, chr_cores = 1, n_bam) {

 chr_idx <- unique(region_df$chrom)

  region_df <-
  mclapply(chr_idx, function(i_chr) {
    res_df <- region_df %>%
      filter(chrom == i_chr) %>%
      mutate(
        id = paste0(chrom, ":", start, "-", end),
        cov = chrCov(f_bam, chrom, start, end, n_bam)) %>%
        separate_rows(cov, sep = ",") %>%
      group_by(id) %>%
      summarise(metric = regionMetric(cov)) %>%
      separate(
        metric, into = c("mean", "median", "rms"),
        sep = ",", convert = TRUE) %>%
      select(id, mean, median, rms)

      message(str_c("Completed run chr:", i_chr))

      return(res_df)
  }, mc.cores = chr_cores) %>%
  bind_rows()

  message(paste("Completed file:", f_bam))
  return(region_df)
}

#' chrIntersect
#'
#' @param region_df
#' @param bam_c
#'
#' @importFrom dplyr filter mutate select if_else
#' @importFrom stringr str_replace
#'
chrIntersect <- function(region_df, bam_c) {
chr_v <- as.character(bam_c)
region_df %>%
  mutate(
    chrom = as.character(chrom),
    chr = if_else(chrom %in% chr_v, chrom, ""),
    chr2 = if_else(paste0("chr", chrom) %in% chr_v,
      paste0("chr", chrom), ""),
    chr3 = if_else(chrom %in% paste0("chr", chr_v), chrom, "")) %>%
  filter(chr != "" | chr2 != "" | chr3 != "") %>%
  mutate(
    chrom = if_else(chr != "", chr, chr2),
    chrom = if_else(chrom == "", str_replace(chr3, "chr", ""), chrom)) %>%
  select(chrom, start, end)
}

chrCov <- function(f_name, chr, start, end, n_read) {

  chr <- unique(chr)

  sig <-
    read_cov(f_name, chr, start, end, n_read)

  sapply(1:length(start), function(i_start) {
    paste0(sig[(start[i_start]):end[i_start]], collapse = ",")
  })

}

#' Calculate metric
#'
#' @param reg_cov
#'
#' @importFrom wavethresh wd threshold
#'
regionMetric <- function(reg_cov) {
  if (anyNA(reg_cov)) {
    "NA,NA,NA"
  } else {
    wd_sig <- wd(reg_cov, filter.number = 1, family = "DaubExPhase")
    wd_sig_tr <- threshold(wd_sig, by.level = TRUE,
        policy = "universal", return.thresh = TRUE)
    paste(spectorRa(wd_sig_tr), spectorMed(wd_sig_tr),
      spectorRms(wd_sig_tr), sep = ",")

  }
}


spectorRa <- function(wd_thr) {
  mean(abs(wd_thr), na.rm = T)
}

spectorRms <- function(wd_thr) {
  sqrt(mean(wd_thr^2))
}

spectorMed <- function(wd_thr) {
  median(abs(wd_thr), na.rm = T)
}

sigNorm <- function(signal) {
  if (max(signal, na.rm = TRUE) != 0) {
    signal <- (signal - min(signal, na.rm = TRUE)) /
              (max(signal, na.rm = TRUE) - min(signal, na.rm = TRUE))
  }

  return(signal)
}
