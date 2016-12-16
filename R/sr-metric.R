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
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#'
#' @export
#'
spector_metric <- function(f_bam = NULL, region_size = NULL, f_bed = NULL,
                          bed.header = FALSE, metric = "wavelet",
                          f.method = NA, region_giab = TRUE) {

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
# EXTRACT BAM FILE STATS USING SAMTOOLS
# ==========================================================================

  bam_stats <- nReadsBam(f_bam)
  n_read <- bam_stats$n_read

# ==========================================================================
# COMPARE CHROMOSOME LISTS BETWEEN BAM AND BED FILE
# ==========================================================================

  chr.i <- .chr_intersect(region_df, bam_stats$chrom)
  message(c("Running spector for:\n", paste("Chrom. ", chr.i, "\n", sep = "")))

#
# Subset region_df by chromosome intersect
# --------------------------------------------------------------------------

  region_df <- region_df %>%
    dplyr::filter(chrom %in% chr.i | chr %in% chr.i) %>%
    # make sure only chr found in bam file are used
    dplyr::filter(!is.na(chrom) & !is.na(chr)) %>%
    dplyr::rowwise()

  if (all(region_df$chrom %in% chr.i)) {
    region_df <-  region_df %>%
      dplyr::do(.region_metric(f_bam, .$chrom, .$start, .$end, n_read,
        metric = metric, methods = f.method)) %>%
      data.frame(chr = region_df$chr, chr.bed = region_df$chrom, id = region_df$id,
        stringsAsFactors = FALSE) %>%
      dplyr::tbl_df()
  } else if (all(region_df$chr %in% chr.i)) {
    region_df <-  region_df %>%
      dplyr::do(.region_metric(f_bam, .$chr, .$start, .$end, n_read,
        metric = metric, methods = f.method)) %>%
      data.frame(chr = region_df$chr, chr.bed = region_df$chr, id = region_df$id,
        stringsAsFactors = FALSE) %>%
      dplyr::tbl_df()
  } else {
    stop("Issue matching chromosome name in bed file with chromosome
      name in bam file")
  }


  if (metric == "fractal") {
    region_df <- region_df %>%
      dplyr::mutate(Df = replace(Df, which(Df == 0), NA))
  } else {
    region_df <- region_df %>%
      dplyr::mutate(R_a = replace(R_a, which(R_a == 0), NA),
              R_rms = replace(R_rms, which(R_rms == 0), NA))
  }


  message(paste("Completed file:", f_bam))
  return(region_df)
}

#
# TODO
#
# => fix check for style of chromosome name

.chr_intersect <- function(region_df, bam.c) {
  if (region_df$chr[1] %in% bam.c) {
    chr <- intersect(region_df$chr, bam.c)
  } else if (is.integer(region_df$chrom)) {
    chr <- intersect(region_df$chrom, bam.c)
  } else if (is.character(region_df$chrom)) {
    chr <- intersect(region_df$chrom, bam.c)
  }else {
    stop("No overlap found for chromosome names between bam and bed files \n
      Check files")
  }

  return(chr)
}

.region_metric <- function(f.name, chr, start, end, n_read, metric = "wavelet",
                            methods = NA, n.lag = "auto",
                            w.size = length(sig)) {

  sig <- read_cov(f.name, chr, start, end, n_read)

  if (metric == "wavelet") {
    #
    # TODO
    #
    # => Make sure this is the best option should all regions be dropped
    #
    if (anyNA(sig)) {
      data.frame(R_a = NA, R_rms = NA)
    } else {
      wd.sig <- wavethresh::wd(sig, filter.number = 1, family = "DaubExPhase")
      wd.sig.tr <- wavethresh::threshold(wd.sig, by.level = TRUE,
          policy = "universal", return.thresh = TRUE)
      data.frame(
        R_a = .spector_ra(wd.sig.tr),
        R_rms = .spector_rms(wd.sig.tr),
        R_md = .spector_med(wd.sig.tr)
        )
    }
  }
}


.spector_ra <- function(wd.thr) {
  mean(abs(wd.thr), na.rm = T)
}

.spector_rms <- function(wd.thr) {
  sqrt(mean(wd.thr^2))
}

.spector_med <- function(wd.thr) {
  sqrt(median(abs(wd.thr), na.rm = T))
}

