#' Calculate SpECtOR metric for one bam file
#'
#' @param f.bam
#' @param r.region
#' @param f.bed
#' @param bed.header
#' @param metric
#'
#' @return Metric for all specified regions
#' @export
#'
spector_metric <- function(f.bam = NULL, stl_cmd = NULL, r.region = '10k',
                          f.bed = NULL, bed.header = FALSE, metric = "wavelet",
                          f.method = NA) {

  ## Load data frame from bed file --------------------------------------------------
  if (r.region == "custom" & is.null(f.bed)) {
    stop(paste("Regions of genome chosen as 'custom' but no bed file provided"))
  } else if (r.region == "custom" || !is.null(f.bed)) {
    bed.d <- spector_bed(file = f.bed, header = bed.header, ucsc.coord = TRUE)
  } else if (r.region == "10k") {
    bed.d <- spector:::giab.10k
  } else if (r.region == "20k") {
    bed.d <- spector:::giab.20k
  }

### Extract bam file stats using samtools ===========================================
  tmp.bam <- .nReadsBam(f.bam, cmd = stl_cmd)
  n.read <- tmp.bam$n.read

### Compare chromosome lists between bam and bed file ===============================
  chr.i <- .chr_intersect(bed.d, tmp.bam$chrom)
  message(c("Running spector for:\n", paste("Chrom. ", chr.i, "\n", sep = "")))

  ## Subset bed.d by chromosome intersect -------------------------------------------
  bed.d <- bed.d %>%
    dplyr::filter(chrom %in% chr.i | chr %in% chr.i) %>%  # make sure only chr found in bam file are used
    dplyr::filter(!is.na(chrom) & !is.na(chr)) %>%
    dplyr::rowwise()

  if (bed.d$chrom %in% chr.i) {
    bed.d <-  bed.d %>%
      dplyr::do(.region_metric(f.bam, .$chrom, .$start, .$end, n.read,
        metric = metric, methods = f.method)) %>%
      data.frame(chr = bed.d$chr, chr.bed = bed.d$chrom, id = bed.d$id) %>%
      dplyr::tbl_df()
  } else if (bed.d$chr %in% chr.i ) {
    bed.d <-  bed.d %>%
      dplyr::do(.region_metric(f.bam, .$chr, .$start, .$end, n.read,
        metric = metric, methods = f.method)) %>%
      data.frame(chr = bed.d$chr, chr.bed = bed.d$chr, id = bed.d$id) %>%
      dplyr::tbl_df()
  } else {
    stop("Issue matching chromosome name in bed file with chromosome name in bam file")
  }


  if (metric == "fractal") {
    bed.d <- bed.d %>%
      dplyr::mutate(Df = replace(Df, which(Df == 0), NA))
  } else {
    bed.d <- bed.d %>%
      dplyr::mutate(R_a = replace(R_a, which(R_a == 0), NA),
              R_rms = replace(R_rms, which(R_rms == 0), NA))
  }

  bed.d$chrom <- as.character(bed.d$chrom)

  message(paste("Completed file:", f.bam))
  return(bed.d)
}

.chr_intersect <- function(bed.d, bam.c) {
  if (bed.d$chr %in% bam.c) {
    chr <- intersect(bed.d$chr, bam.c)
  } else if (is.integer(bed.d$chrom)) {
    chr <- intersect(bed.d$chrom, bam.c)
  } else if (is.character(bed.d$chrom)) {
    chr <- intersect(bed.d$chrom, bam.c)
  }else {
    stop("No overlap found for chromosome names between bam and bed files \nCheck files")
  }

  return(chr)
}

.region_metric <- function(f.name, chr, start, end, n.read, metric = "wavelet",
                            methods = NA, n.lag = "auto",
                            w.size = length(sig)) {

  sig <- read_cov(f.name, chr, start, end, n.read)

  if (metric == "wavelet") {
    if (is.na(sig[1])) {
      data.frame(R_a = NA, R_rms = NA)
    } else {
      wd.sig <- wavethresh::wd(sig, filter.number = 1, family = "DaubExPhase")
      wd.sig.tr <- wavethresh::threshold(wd.sig, by.level = TRUE,
          policy = "universal", return.thresh = TRUE)
      data.frame(R_a = .spector_ra(wd.sig.tr), R_rms = .spector_rms(wd.sig.tr))
    }
  }
}


.spector_ra <- function(wd.thr) {
  mean(abs(wd.thr))
}

.spector_rms <- function(wd.thr) {
  sqrt(mean(wd.thr^2))
}
