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
spector_metric <- function(f_bam = NULL, stl_cmd = NULL, region_size = NULL,
                          f_bed = NULL, bed.header = FALSE, metric = "wavelet",
                          f.method = NA, region_giab = TRUE) {

#
# Load data frame from bed file
# --------------------------------------------------------------------------
  if (region_giab) {
    bed.d <-
      giab_10k %>%
        mutate(reg_length = end - start) %>%
        bed_region_split(region_size)

  } else if (!region_giab) {
    bed.d <- read_bed(bed_file = f_bed, header = bed.header,
      region_size = region_size)
  }

# ==========================================================================
# EXTRACT BAM FILE STATS USING SAMTOOLS
# ==========================================================================

  tmp.bam <- .nReadsBam(f_bam, cmd = stl_cmd)
  n.read <- tmp.bam$n.read

# ==========================================================================
# COMPARE CHROMOSOME LISTS BETWEEN BAM AND BED FILE
# ==========================================================================

  chr.i <- .chr_intersect(bed.d, tmp.bam$chrom)
  message(c("Running spector for:\n", paste("Chrom. ", chr.i, "\n", sep = "")))

#
# Subset bed.d by chromosome intersect
# --------------------------------------------------------------------------

  bed.d <- bed.d %>%
    dplyr::filter(chrom %in% chr.i | chr %in% chr.i) %>%
    # make sure only chr found in bam file are used
    dplyr::filter(!is.na(chrom) & !is.na(chr)) %>%
    dplyr::rowwise()

  if (all(bed.d$chrom %in% chr.i)) {
    bed.d <-  bed.d %>%
      dplyr::do(.region_metric(f_bam, .$chrom, .$start, .$end, n.read,
        metric = metric, methods = f.method)) %>%
      data.frame(chr = bed.d$chr, chr.bed = bed.d$chrom, id = bed.d$id,
        stringsAsFactors = FALSE) %>%
      dplyr::tbl_df()
  } else if (all(bed.d$chr %in% chr.i)) {
    bed.d <-  bed.d %>%
      dplyr::do(.region_metric(f_bam, .$chr, .$start, .$end, n.read,
        metric = metric, methods = f.method)) %>%
      data.frame(chr = bed.d$chr, chr.bed = bed.d$chr, id = bed.d$id,
        stringsAsFactors = FALSE) %>%
      dplyr::tbl_df()
  } else {
    stop("Issue matching chromosome name in bed file with chromosome
      name in bam file")
  }


  if (metric == "fractal") {
    bed.d <- bed.d %>%
      dplyr::mutate(Df = replace(Df, which(Df == 0), NA))
  } else {
    bed.d <- bed.d %>%
      dplyr::mutate(R_a = replace(R_a, which(R_a == 0), NA),
              R_rms = replace(R_rms, which(R_rms == 0), NA))
  }


  message(paste("Completed file:", f_bam))
  return(bed.d)
}

#
# TODO
#
# => fix check for style of chromosome name

.chr_intersect <- function(bed.d, bam.c) {
  if (bed.d$chr[1] %in% bam.c) {
    chr <- intersect(bed.d$chr, bam.c)
  } else if (is.integer(bed.d$chrom)) {
    chr <- intersect(bed.d$chrom, bam.c)
  } else if (is.character(bed.d$chrom)) {
    chr <- intersect(bed.d$chrom, bam.c)
  }else {
    stop("No overlap found for chromosome names between bam and bed files \n
      Check files")
  }

  return(chr)
}

.region_metric <- function(f.name, chr, start, end, n.read, metric = "wavelet",
                            methods = NA, n.lag = "auto",
                            w.size = length(sig)) {

  sig <- read_cov(f.name, chr, start, end, n.read)

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

