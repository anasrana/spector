#' @importFrom dplyr mutate group_by summarise bind_rows filter is.tbl
#' @importFrom stringr str_c
#' @importFrom magrittr %>%
#' @importFrom tidyr separate
#' @importFrom parallel mclapply
#' @importFrom purrr map map_dbl
#'
spectorMetric <- function(region_df, f_bam = NULL, chr_cores = 1, n_bam,
                          met = smr) {

if (!is.tbl(region_df)) {
  stop("region_df needs to be a tbl_df object after loading\n",
    "Something went wrong", call. = FALSE)
}

 chr_idx <- unique(region_df$chrom)

  region_df <-
  mclapply(chr_idx, function(i_chr) {

    res_df <- region_df %>%
      regionCovDf(chr = i_chr, bam_file = f_bam, bamReadCount = n_bam) %>%
      mutate(wd_thresh = map(cov, regionMetric),
             metric = map_dbl(wd_thresh, spectorRms),
             region.status = if_else(is.infinite(metric),
                                     "not.covered", "covered"),
             metric = if_else(is.infinite(metric), as.double(NA), metric))

  if (met == "all") {
    res_df <-
    res_df %>%
      mutate(median = map_dbl(wd_thresh, spectorMed),
             mean = map_dbl(wd_thresh, spectorRa)) %>%
      select(chrom, start, end, metric, mean, median, region.status)
    } else if (met == "rms") {
      res_df <-
      res_df %>%
        select(chrom, start, end, metric, region.status)
    }

      message(str_c("Completed run chr:", i_chr))

      return(res_df)
  }, mc.cores = chr_cores) %>%
  bind_rows()

  message(paste("Completed file:", f_bam))
  return(region_df)
}

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

  lapply(seq_along(start), function(i_start) {
    sig[(start[i_start]):end[i_start]]
  })

}

#' @importFrom wavethresh wd threshold
#'
regionMetricAll <- function(reg_cov) {
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

#' @importFrom wavethresh wd threshold
#'
regionMetric <- function(reg_cov) {

  if (anyNA(reg_cov)) {
    NA
  } else {
    wd_sig <- wd(as.vector(reg_cov), filter.number = 1, family = "DaubExPhase")
    wd_sig_tr <- threshold(wd_sig, by.level = TRUE,
        policy = "universal", return.thresh = TRUE)

    return(wd_sig_tr)
  }
}

spectorRa <- function(wd_thr) {
  1 / mean(abs(wd_thr), na.rm = T)
}

spectorRms <- function(wd_thr) {
  1 / sqrt(mean(wd_thr^2))
}

spectorMed <- function(wd_thr) {
  1 / median(abs(wd_thr), na.rm = T)
}

sigNorm <- function(signal) {
  if (max(signal, na.rm = TRUE) != 0) {
    signal <- (signal - min(signal, na.rm = TRUE)) /
              (max(signal, na.rm = TRUE) - min(signal, na.rm = TRUE))
  }

  return(signal)
}

#' @importFrom tidyr separate_rows
#' @importFrom dplyr mutate filter
#'
regionCovDf <- function(region_df, chr, bam_file, bamReadCount) {
  region_df %>%
      filter(chrom == chr) %>%
      mutate(
        id = paste0(chrom, ":", start, "-", end),
        cov = chrCov(bam_file, chrom, start, end, bamReadCount))
}
