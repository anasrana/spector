#' Compute QC (spector metric) for bam files.
#'
#' Wavelet based technique to compute a quality metric for regions
#' across the genome. The `spector_qc()` is the recommended function to
#' access `spector` since it includes several checks and allows for more
#' flexibility than some of the downstream functions.
#'
#' The `spector_qc` function is the main function to use for QC in the
#' `spector` package. It will compute a quality control metric for specific
#' regions across the genome.
#' The default regions, supplied in the package, are based
#' on the genome in a bottle project
#' ([giab](http://jimb.stanford.edu/giab/)) reliable regions, calculated
#' using `ReliableGenome` ([RG](http://github.com/popitsch/wtchg-rg)).
#' It is also possible to supply custom regions as a bed file or a
#' `data.frame` object.
#'
#' It is important to supply full paths to `f_bam`, and `f_bed`.
#' Though the path can be relative to the current working directory, which can
#' be set with [base::setwd()]. This also applies to the first
#' column of a parameter file that can be supplied to `f_bam`.
#'
#' @param f_bam a string with path to `bam` file(s), it can link to
#'        `*.bam` file (the full relative path is required), a folder with
#'        `*.bam` files, or a file with with structure specified later.
#' @param f_bed file path for bed file to override default giab.
#' @param region_giab logical. Indicates whether or not giab regions are used,
#'        defaults to `TRUE`.
#' @param region_size integer. Choose size of regions to calculate metric value.
#'        The default `= NULL` means `region_size = ` maximum power
#'        of 2 that fits in the smallest region.
#' @param file_type type of file passed to f_bam (Optional). This is to ensure
#'        the automated checks pick up the correct format. The possible options
#'        are `"list"`, `"bam"`, and `"dir"`.
#' @param out_F (Optional) Folder path to save output. If omitted results
#'        will be returned, but not saved.
#' @param save_out logical. Indicating if output from `spector_qc()`
#'        should be saved.
#' @param silent logical. Default `FALSE`, if `TRUE` there is no
#'        progress update for the code.
#' @param smr_var deprecated. Variable to use to compute summary.
#' @param file_cores integer. Optional number indicating if the QC should be
#'        computed in parallel across all input files.
#' @param bed_header logical. `TRUE` if bed file has a header, the default
#'        value is `FALSE`.
#' @param chr_cores integer. Optional number indicating if the QC should be
#'        computed in parallel across chromosomes. Default value is `1`.
#'
#' @return Output is a `tbl_df` object with a metric value for each
#'         region. Optionally the output can also be saved to file, but only if
#'         `out_F` is provided.
#'
#' @examples
#' \dontrun{
#' # Compute QC on sampl1.bam with default options
#' spector_qc(f_bam = "sample1.bam")
#'}
#'
#' s1_path <- spector_sample("sample1.bam")
#' basic_path <- spector_sample("basic.bed")
#'
#' # Compute QC on sample1.bam with custom region size
#' spector_qc(f_bam = s1_path, region_size = 2^14)
#'
#' # Compute QC on sample1.bam with custom bed file
#' spector_qc(f_bam = s1_path, f_bed = basic_path)
#'
#' # Compute QC and save output results
#' spector_qc(f_bam = s1_path, f_bed = basic_path, out_F = "~/",
#'  save_out = TRUE)
#'
#' @export
#'
#' @importFrom stringr str_c
#' @importFrom utils capture.output count.fields
#'
spector_qc <- function(f_bam = NULL, f_bed = NULL, region_giab = TRUE,
                       region_size = NULL, file_type = "bam", out_F = NULL,
                       save_out = FALSE, silent = FALSE, smr_var = "rms",
                       file_cores = 1, chr_cores = 1, bed_header = FALSE) {

# ==========================================================================
# CHECKING VARIABLE
# ==========================================================================

## Check if have read access to f_bam
  if (file.access(f_bam, mode = 4) == -1) {
    stop(str_c("'", f_bam, "'\nYou do not have access to the file or the",
               "file does not exist"), call. = FALSE)
  }

## Create output folder if it doesn't exist already
  if (save_out) {

    out_F <- outTrailingFix(out_F)

    if (!dir.exists(out_F)) {
      dir.create(out_F)

      message(str_c("Folder '", out_F,
                    "' created in current working directory ('", getwd(), "')"))
    }
  }

## should output be printed

  if (silent) {
    output_capture = "message"
  } else {
    output_capture = NULL
  }

## clash between parallel

if (file_cores > 1 & chr_cores > 1) {
  stop("It is advisable to only run in parallel across chromosomes OR files!",
    "\n\tThey are nested functions.", call. = FALSE)
}

# ==========================================================================
# Computing metric after checking what f_bam is
# ==========================================================================

  if (dir.exists(f_bam)) {
    fs_bam <- paste(f_bam, '/', list.files(path = f_bam, pattern = "*.bam$"),
      sep = "")

    id_bam <- gsub(".bam","", x = basename(fs_bam))

    capture.output(
      srm_df <- spectorList(fs_bam, id_v = id_bam, s_v = NULL,
                            out_F = out_F, file_cores = file_cores,
                            chr_cores = chr_cores, save_out = save_out,
                            region_giab = region_giab,
                            region_size = region_size, f_bed = f_bed,
                            bed_header = bed_header),
      type = output_capture)

  } else if (file.exists(f_bam)) {
    if (file_type == "list") {

      bam_pars <- read_par_file(id_path = f_bam)
      unpackList(bam_pars)

      capture.output(
        srm_df <- spectorList(fs_bam = fs_bam, id_v = id_bam, s_v = sample_type,
                              out_F = out_F, file_cores = file_cores,
                              chr_cores = chr_cores, save_out = save_out,
                              region_giab = region_giab,
                              region_size = region_size, f_bed = f_bed,
                              bed_header = bed_header),
        type = output_capture)


    } else if (grep("*.bam$", x = f_bam) > 0 | file_type == "bam") {

      capture.output(
        srm_df <- spectorFile(f_bam, out_F = out_F, save_out = save_out,
                              chr_cores = chr_cores, region_giab = region_giab,
                              region_size = region_size, f_bed = f_bed,
                              header = bed_header),
        type = output_capture)

    } else if(!file.exists(f_bam)) {

      stop(str_c("f_bam: ", f_bam, " - file/folder not found",
                 "\nMake sure you provide the correct path"), call. = FALSE)
    }
  }

# ==========================================================================
# Save output and return
# ==========================================================================

  if (save_out) {
    # save outputs
    saveMerged(res_v = srm_df, out = out_F)

    saveSummary(res = srm_df, out = out_F, var_s = smr_var)
  }

  return(srm_df = srm_df)
}

#
# Auxillary function
# --------------------------------------------------------------------------

unpackList <- function(object) {
  for(.x in names(object)){
    assign(value = object[[.x]], x=.x, envir = parent.frame())
  }
}

#' @importFrom dplyr group_by summarise
#' @importFrom utils write.csv
#'
saveSummary <- function(res, var_s, out) {

  if (var_s == "rms") {
    stat_spector <- res %>%
      group_by(id_bam) %>%
      summarise(mean_rm = mean(1 / rms, na.rm = TRUE),
                       median_rm = median(1 / rms, na.rm = TRUE),
                       sd_rm = sd(1 / rms, na.rm = TRUE),
                       iqr_rm = IQR(1 / rms, na.rm = TRUE))
  } else if (var_s == "mean") {
    stat_spector <- res %>%
      group_by(id) %>%
      summarise(mean_rm = mean(1 / mean, na.rm = TRUE),
                       median_rm = median(1 / mean, na.rm = TRUE),
                       sd_rm = sd(1 / mean, na.rm = TRUE),
                       iqr_rm = IQR(1 / mean, na.rm = TRUE))
  }

  write.csv(stat_spector,
            file = paste(out, "SUMMARY_STAT_metric.csv", sep = ""),
            row.names = FALSE)
}

#' @importFrom utils write.csv
#'
saveMerged <- function(res_v, out) {
    out_path <- paste(out, "results_bam_out_merged.csv", sep = "")
    write.csv(res_v, file = out_path, row.names = FALSE)
}

#' Read spector id file.
#'
#' @param id_path path to the id file
#'
#' @export
#'
#' @importFrom tibble as_data_frame
#' @importFrom utils read.table count.fields
#' @family data import functions
#'
read_par_file <- function(id_path) {
  list_names <- c("fs_bam", "id_bam", "sample_type")

  no_cols <- max(count.fields(id_path, sep = "\t"))

  id_df <-   read.table(id_path, strip.white = TRUE, stringsAsFactors = FALSE,
    col.names = c(list_names, rep(NA, no_cols - length(list_names))),
    colClasses = c(rep("character", length(list_names)),
      rep("NULL", no_cols - length(list_names)))) %>%
    tibble::as_data_frame()

  return(id_df)
}

outTrailingFix <- function(out_F) {
  tmp_out <- unlist(strsplit(out_F, split = ""))
  if (!(tmp_out[length(tmp_out)] == "/")) {
    out_F <- paste(out_F, "/", sep = "")
  }
  return(out_F)

}

#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
