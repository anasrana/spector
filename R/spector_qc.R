#' Compute QC (spector metric) for bam files.
#'
#' Wavelet based technique to compute a quality metric for regions
#' across the genome.
#'
#' This is the main function to use for QC in the \pkg{spector} package. It will
#' compute a quality control metric for specific regions across the genome.
#' The default regions, supplied in the package, are based
#' on the genome in a bottle project
#' (\href{http://jimb.stanford.edu/giab/}{giab}) reliable regions, calculated
#' using \code{ReliableGenome} (\href{http://github.com/popitsch/wtchg-rg}{RG}).
#' It is also possible to supply custom regions as a bed file or a
#' \code{data.frame} object.
#'
#' It is important to supply full paths to \code{f_bam}, and \code{f_bed}.
#' Though the path can be relative to the current working directory, which can
#' be set with \code{\link[base]{setwd}()}. This also applies to the first
#' column of a parameter file that can be supplied to \code{f_bam}.
#'
#' @param f_bam a string with path to \code{bam} file(s), it can link to
#'        \code{*.bam} file (the full relative path is required), a folder with
#'        \code{*.bam} files, or a file with with structure specified later.
#' @param file_type type of file passed to f_bam (Optional). This is to ensure
#'        the automated checks pick up the correct format. The possible options
#'        are \code{"list"}, \code{"bam"}, and \code{"dir"}.
#' @param f_delim the delimiter character used in \code{f_bam} file if not
#'        \code{bam} or folder (Optional). It is only used when
#'        \code{file_type = "list"}, with the default \code{f_delim = "\t"}.
#' @param out_F (Optional) Folder path to save output. If omitted results
#'        will be returned, but not saved.
#' @param file_cores integer. Optional number indicating if the QC should be
#'        computed in parallel across all input files.
#' @param smr_var deprecated. Variable to use to compute summary.
#' @param save_out logical. Indicating if output from \code{spector_qc()}
#'        should be saved.
#' @param ... variables passed on to downstream function, please see
#'        \code{\link{spector_metric}()} for further details on the parameters.
#'
#' @return Output is a \code{tbl_df} object with a metric value for each region.
#'         Optionally the output can also be saved, if \code{out_F} is provided.
#'
#' @examples
#' \dontrun{
#' # Compute QC on sampl1.bam with default options
#' spector_qc(f_bam = "sample1.bam")
#'}
#'
#' # Compute QC on sample1.bam with custom region size
#' # Note. The full command is bed_region_size, but because of partial matching
#' # we can use this partial command.
#' spector_qc(f_bam = "tests/testthat/sample1.bam", region_size = 2^14)
#'
#' # Compute QC on sample1.bam with custom bed file
#' spector_qc(f_bam = "tests/testthat/sample1.bam",
#'  f_bed = "tests/testthat/basic.bed")
#'
#' # Compute QC and save output results
#' spector_qc(f_bam = "tests/testthat/sample1.bam",
#'  f_bed = "tests/testthat/basic.bed", out_F = "~/", save_out = T)
#'
#' @export
#'
#' @importFrom stringr str_c
#'
spector_qc <- function(f_bam = NULL,
                    file_type = "list",
                    f_delim = "\t",
                    out_F = NULL,
                    file_cores = 1,
                    smr_var = "rms",
                    save_out = FALSE,
                    ...) {

  ## Check if have read access to f_bam
  if (file.access(f_bam, mode = 4) == -1) {
    stop(str_c("'", f_bam,
      "'\nYou do not have acces to the file or the file does not exist"))
  }

  ## Create output folder if it doesn't exist already
  if (save_out) {
    out_F <- outTrailingFix(out_F)
    if (!dir.exists(out_F)) {
      dir.create(out_F)
      message(paste("Folder '", out_F,
        "' created in current working directory ('", getwd(), "')", sep = ""))
    }
  }

#
# Computing metric after checking what f_bam is
# --------------------------------------------------------------------------

  if (dir.exists(f_bam)) {
    fs_bam <- paste(f_bam, '/', list.files(path = f_bam, pattern = "*.bam$"),
      sep = "")
    id_bam <- gsub(".bam","", x = basename(fs_bam))
    srm_df <- spectorList(fs_bam, id_v = id_bam, s_v = NULL,
                           out_F = out_F, file_cores = file_cores, save_out = save_out,
                           ...)

  } else if (file.exists(f_bam)) {
    if (file_type == "list") {

      fs_bam <- readIdFile(id_path = f_bam)
      unpackList(fs_bam)
      srm_df <- spectorList(fs_bam = fs_bam, id_v = id_bam, s_v = sample_type,
                            out_F = out_F, file_cores = file_cores,
                            save_out = save_out, ...)


    } else if (grep("*.bam$", x = f_bam) > 0 | file_type == "bam") {

      srm_df <- spectorFile(f_bam, out_F = out_F, save_out = save_out, ...)

    } else if(!file.exists(f_bam)) {
      stop("f_bam: file / folder not found
            Make sure you provide the correct path")
    }
  }

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

#' saveSummary
#'
#' @param res
#' @param var_s
#' @param out
#'
#' @importFrom dplyr group_by summarise
#'
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

saveMerged <- function(res_v, out) {
    out_path <- paste(out, "results_bam_out_merged.csv", sep = "")
    write.csv(res_v, file = out_path, row.names = FALSE)
}

#' Read spector id file.
#'
#' @param id_path path to the id file
#'
#' @importFrom tibble as_data_frame
#' @importFrom utils read.table
#'
readIdFile <- function(id_path) {
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
