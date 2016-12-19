#' Calculate SpECtOR metric for one bam file
#'
#' @param bed.header binary, \code{TRUE} if bed file has a header
#' @param f_bam
#' @param f_bed file path for bed file if \code{region_giab = FALSE}
#' @param region_giab logial indicating if giab regions use or not, the default
#'                    TRUE
#' @param region_size choose size of the regions defualt NULL chooses the max
#'    power of 2 that fits in the smallest region
#'
#' @return Metric for all specified regions
#' @importFrom dplyr mutate group_by summarise
#' @importFrom magrittr %>%
#' @importFrom tidyr separate separate_rows
#'
#' @export
#'
spector_metric <- function(f_bam = NULL, region_size = NULL, f_bed = NULL,
                          bed.header = FALSE, f.method = NA,
                          region_giab = TRUE) {

#
# Load data frame from bed file
# --------------------------------------------------------------------------
  if (region_giab) {
    region_df <-
      giab_10k %>%
        mutate(reg_length = end - start) %>%
        bed_region_split(region_size)

  } else if (!region_giab) {
    region_df <- read_bed(bed_file = f_bed, header = bed.header,
      region_size = region_size)
  }

# ==========================================================================
# EXTRACT BAM FILE STATS USING RSamtools
# ==========================================================================

  bam_stats <- nReadsBam(f_bam)
  n_read <- bam_stats$n_read

# ==========================================================================
# COMPARE CHROMOSOME LISTS BETWEEN BAM AND BED FILE
# ==========================================================================

  region_df <- chrIntersect(region_df, bam_stats$chrom)

  if (nrow(region_df) > 0) {
    message(c("Running spector for:\n",
      paste0("Chromosome: ", unique(region_df$chrom), "\n")))
  } else {
    stop("Issues matching chr in '*.bed' and '*.bam' or no overlap")
  }

#
# Subset region_df by chromosome intersect
# --------------------------------------------------------------------------

  region_df <-
    region_df %>%
    group_by(chrom) %>%
    mutate(
      id = paste0(chrom, ":", start, "-", end),
      cov = get_chr_cov(f_bam, chrom, start, end, n_read)) %>%
    separate_rows(cov, sep = ",") %>%
    group_by(id) %>%
    summarise(metric = region_metric(cov)) %>%
    separate(
      metric, into = c("mean", "median", "rms"),
      sep = ",", convert = TRUE) %>%
    select(id, mean, median, rms)

  message(paste("Completed file:", f_bam))
  return(region_df)
}

#' chrIntersect
#'
#' @param region_df
#' @param bam_c
#'
#' @importFrom dplyr filter mutate select if_else
#'
chrIntersect <- function(region_df, bam_c) {
chr_v <- as.character(bam_c)
region_df %>%
  mutate(
    chr = if_else(chrom %in% chr_v, chrom, ""),
    chr2 = if_else(paste0("chr", chrom) %in% chr_v,
      paste0("chr", chrom), "")) %>%
  filter(chr != "" | chr2 != "") %>%
  mutate(chrom = if_else(chr != "", chr, chr2)) %>%
  select(chrom, start, end)
}

get_chr_cov <- function(f.name, chr, start, end, n_reawd.sig.trd) {

  chr <- unique(chr)

  sig <-
    read_cov(f.name, chr, start, end, n_read)

  sapply(1:length(start), function(i_start) {
    paste0(sig[(start[i_start] + 1):end[i_start]], collapse = ",")
  })

}

region_metric <- function(reg_cov) {
  if (anyNA(reg_cov)) {
    "NA,NA,NA"
  } else {
    wd.sig <- wavethresh::wd(reg_cov, filter.number = 1, family = "DaubExPhase")
    wd.sig.tr <- wavethresh::threshold(wd.sig, by.level = TRUE,
        policy = "universal", return.thresh = TRUE)
    paste(.spector_ra(wd.sig.tr), .spector_med(wd.sig.tr),
      .spector_rms(wd.sig.tr), sep = ",")

  }

}


.spector_ra <- function(wd.thr) {
  mean(abs(wd.thr), na.rm = T)
}

.spector_rms <- function(wd.thr) {
  sqrt(mean(wd.thr^2))
}

.spector_med <- function(wd.thr) {
  median(abs(wd.thr), na.rm = T)
}

.sig_norm <- function(signal) {
  if (max(signal, na.rm = TRUE) != 0) {
    signal <- (signal - min(signal, na.rm = TRUE)) /
              (max(signal, na.rm = TRUE) - min(signal, na.rm = TRUE))
  }

  return(signal)
}
