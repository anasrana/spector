#' @importFrom dplyr is.tbl
#' @importFrom stringr str_c
#'
spectorFile <- function(f_bam, id_bam = NULL, s_prep = NULL, out_F = NULL,
                        save_out, chr_cores, region_giab = TRUE,
                        region_size = NULL, f_bed = NULL,
                        header = FALSE, smr = "rms") {

  message(paste("Running on file:", f_bam, "\n=>\n"))

  #
  # Load region data frame from bed file or filter giab regions
  # --------------------------------------------------------------------------

  if(!is.tbl(f_bed)) {

    region_df <- getRegions(region_giab = region_giab, f_bed = f_bed,
                            region_size = region_size, header = header)
  } else {
    region_df <- f_bed
  }


  # extract bam file stats using rsamtools
  # --------------------------------------------------------------------------

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
    stop("Issues matching chr in '*.bed' and '*.bam' or no overlap",
         call. = FALSE)
  }

  message(str_c(format(nrow(region_df), big.mark = ","),
    " regions identified."))


  srm_df <- spectorMetric(region_df = region_df, f_bam = f_bam,
                           chr_cores = chr_cores, n_bam = n_read, met = smr)
  if (is.null(id_bam)) {
    id_bam <- gsub(".bam", "", x = basename(f_bam))
  }

  srm_df$id_bam <- id_bam
  srm_df$prep <- s_prep


  if (is.null(out_F) & save_out) {
    warning("Output folder not specified, output will not be saved",
      call. = FALSE, domain = 'spector_run')
  } else if (save_out | !is.null(out_F)) {
    out_path <- paste(out_F, id_bam, "_out.csv", sep = "")
    write.csv(srm_df, file = out_path, row.names = FALSE)
    message(paste("Completed, output written to:", out_path))
  }

  return(srm_df)
}

#' @importFrom parallel mclapply
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_c
#'
spectorList <- function(fs_bam, id_v, s_v, out_F, file_cores = 1,
                        save_out, chr_cores, region_giab = TRUE,
                        region_size = NULL, f_bed = NULL,
                        bed_header = FALSE, smr = "rms") {

  # function to be run inside of mclappy
  f_idx <- 1:length(fs_bam)

  srm_df <-
  mclapply(X = f_idx, function(idx) {
    message(str_c("Running file: ", fs_bam[idx]))
      spectorFile(f_bam = fs_bam[idx], id_bam = id_v[idx],
                  s_prep = s_v[idx], out_F = out_F, save_out = save_out,
                  chr_cores = chr_cores,  region_giab = region_giab,
                        region_size = region_size, f_bed = f_bed,
                        header = bed_header, smr = smr)},
      mc.cores = file_cores) %>%
    bind_rows()

  return(srm_df)
}

