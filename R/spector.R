#' Run spector for a specified path or specified files
#'
#' Takes in a path to a folder of .bam files, a single .bam file or a file
#' containing a list with the first column being a file paths and the
#' second column an id and the final column a groupID
#'
#' @param bam_f path to bam file, txt file with paths or folder
#' @param file_type type of file passed to bam_f, options are \code{"list"},
#'                  \code{"bam"}, and \code{"dir"}
#' @param f_delim file delimter, only used when \code{file_type = "list"}
#' @param f_head binary does the \code{bam_f} file have a header,
#'               only used when \code{file_type = "list"}
#' @param out_F path to output folder
#' @param n_core number of cores that should be used
#' @inheritParams spector_metric
#'
#' @return Output files saved in folder \code{out_F}, also saves a reference
#'          file with file names
#' @export
#'
#' @importFrom stringr str_c
#'
spector <- function(bam_f = NULL,
                    file_type = "list",
                    f_delim = "\t",
                    f_head = FALSE,
                    out_F = NULL,
                    n_core = NULL,
                    ...) {

  ## Check if have read access to bam_f
  if (file.access(bam_f, mode = 4) == -1) {
    stop(str_c("'", bam_f,
      "'\nYou do not have acces to the file or the file does not exist"))
  }

  ## Create output folder if it doesn't exist already
  out_F <- outTrailingFix(out_F)
  if (!dir.exists(out_F)) {
    dir.create(out_F)
    message(paste("Folder '", out_F,
      "' created in current working directory ('", getwd(), "')", sep = ""))
  }

#
# Computing metric after checking what bam_f is
# --------------------------------------------------------------------------

  if (dir.exists(bam_f)) {
    fs_bam <- paste(bam_f, '/', list.files(path = bam_f, pattern = "*.bam$"),
      sep = "")
    id_bam <- gsub(".bam","", x = basename(fs_bam))
    srm_df <- spectorList(fs_bam, id_v = id_bam, grp_v = NULL, s_v = NULL,
                           out_F = out_F, n_core = n_core, ...)
    # save outputs
    saveMerged(res_v = srm_df, out = out_F)

  } else if (file.exists(bam_f)) {
    if (file_type == "list") {

      fs_bam <- readIdAssign(id_path = bam_f, f_head = f_head)
      unpackList(fs_bam)
      srm_df <- spectorList(fs_bam = fs_bam, id_v = id_bam, grp_v = gr_bam,
                              out_F = out_F, n_core = n_core, ...)

      # save output
      saveMerged(res_v = srm_df, out = out_F)

    } else if (grep("*.bam$", x = bam_f) > 0 | file_type == "bam") {

      srm_df <- spectorFile(bam_f, out_F = out_F,  ...)

    } else if(!file.exists(bam_f)) {
      stop("bam_f: file / folder not found
            Make sure you provide the correct path")
    }
  }

  saveSummary(res = srm_df, out = out_F)

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

saveSummary <- function(res, var_s, out) {

  if (var_s == "rms") {
    stat_spector <- res %>%
      dplyr::group_by(id_bam) %>%
      dplyr::summarise(mean_rm = mean(1 / rms, na.rm = TRUE),
                       median_rm = median(1 / rms, na.rm = TRUE),
                       sd_rm = sd(1 / rms, na.rm = TRUE),
                       iqr_rm = IQR(1 / rms, na.rm = TRUE))
  } else if (var_s == "mean") {
    stat_spector <- res %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(mean_rm = mean(1 / mean, na.rm = TRUE),
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

readIdAssign <- function(id_path, f_head = FALSE) {
  list.name <- c("fs_bam", 'id_bam', 'gr_bam', 'baseline', 's_prep')
  for (i in 1:length(list.name)) {
    assign(list.name[i], NULL)
  }
  id_df <- dplyr::tbl_df(read.csv(id_path, sep = ",", header = f_head,
    stringsAsFactors = FALSE, strip.white = TRUE))
  colnames(id_df) <- list.name[1:ncol(id_df)]
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
