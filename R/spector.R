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
#' @param plot_type Vector for the type of plots to be saved. Options include
#'                  \code{"boxplot", "circle", "dot"}, see \code{spector_plot} for
#'                  further details.
#' @param plot_var The variable to be used for plotting, default is \code{"rms"}
#'                  for wavelets.
#' @param samtools_cmd System command for samtools including full path if necessary
#' @param r.region choose between the default \code{"10k"} **or** \code{"20k"}
#'                  \code{"custom"} regions
#' @param f.bed file path for bed file if \code{r.region = "custom"}
#' @param bed.header binary, \code{TRUE} if bed file has a header
#' @param ... additional options passed to \code{gr_metric()}

#'
#' @return Output files saved in folder \code{out_F}, also saves a reference file
#'          with file names
#' @export
#'
spector <- function(bam_f = NULL,
                    file_type = "list",
                    f_delim = "\t",
                    f_head = FALSE,
                    out_F = NULL,
                    n_core = NULL,
                    spector_plot = TRUE,
                    plot_type = "boxplot",
                    plot_var = "rms",
                    samtools_cmd = "samtools",
                    ...) {
  ## Create output folder if it doesn't exist already
  out_F <- .out_trailing_fix(out_F)
  if (!dir.exists(out_F)) {
    dir.create(out_F)
    message(paste("Folder '", out_F, "' created in current working directory ('",
                   getwd(), "')", sep = ""))
  }

# ================================================================================
# Computing metric after checking what bam_f is
  if (dir.exists(bam_f)) {
    # TODO include path check
    fs_bam <- paste(bam_f, '/', list.files(path = bam_f, pattern = "*.bam$"),
      sep = "")
    id_bam <- gsub(".bam","", x = basename(fs_bam))
    srm.df <- .spector_list(fs_bam, id_v = id_bam, grp_v = NULL, s_v = NULL,
                           out_F = out_F, n_core = n_core,
                           stl_cmd = samtools_cmd, ...)
    # save outputs
    .save_merged(res_v = srm.df, out = out_F)
    spector_pl <- spector_plot(id_file = NULL, res_df = srm.df, res_p = NULL,
                              f_head = f_head, plot_type = plot_type,
                              plot_var =  plot_var, out_F = out_F)

  } else if (file.exists(bam_f)) {
    if (file_type == "list") {

      fs_bam <- .read_id_assign(id_path = bam_f, f_head = f_head)
      .unpack_list(fs_bam)
      srm.df <- .spector_list(fs_bam = fs_bam, id_v = id_bam, grp_v = gr_bam,
                              out_F = out_F, n_core = n_core,
                              stl_cmd = samtools_cmd, ...)

      # save output
      .save_merged(res_v = srm.df, out = out_F)
      spector_pl <- spector_plot(id_file = bam_f, res_df = srm.df, res_p = NULL,
                                  f_head = f_head, plot_type = plot_type,
                                  plot_var =  plot_var, out_F = out_F)

    } else if (grep("*.bam$", x = bam_f) > 0 | file_type == "bam") {
      srm.df <- .spector_file(bam_f, out_F = out_F, stl_cmd = samtools_cmd, ...)
      # spector_pl <- spector_plot(id_file = NULL, res_df = srm.df, res_p = NULL,
      #                            f_head = FALSE, plot_type = plot_type,
      #                            plot_var =  plot_var, out_F = out_F)
    } else if(!file.exists(bam_f)) {
      stop("bam_f: file / folder not found
            Make sure you provide the correct path")
    }
  }

  .save_summary(res = srm.df, var_s = plot_var, out = out_F)

  return(srm.df = srm.df)
}

#=================================================================================
# Auxillary function


.unpack_list <- function(object) {
  for(.x in names(object)){
    assign(value = object[[.x]], x=.x, envir = parent.frame())
  }
}

.save_summary <- function(res, var_s, out) {

  if (var_s == "rms") {
    stat_spector <- res %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(mean_rm = mean(1 / R_rms, na.rm = TRUE),
                       median_rm = median(1 / R_rms, na.rm = TRUE),
                       sd_rm = sd(1 / R_rms, na.rm = TRUE),
                       iqr_rm = IQR(1 / R_rms, na.rm = TRUE))
  } else if (var_s == "mean") {
    stat_spector <- res %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(mean_rm = mean(1 / R_a, na.rm = TRUE),
                       median_rm = median(1 / R_a, na.rm = TRUE),
                       sd_rm = sd(1 / R_a, na.rm = TRUE),
                       iqr_rm = IQR(1 / R_a, na.rm = TRUE))
  } else if (var_s == "df") {
    stat_spector <- res %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(mean_rm = mean(Df, na.rm = TRUE),
                       median_rm = median(Df, na.rm = TRUE),
                       sd_rm = sd(Df, na.rm = TRUE),
                       iqr_rm = IQR(Df, na.rm = TRUE))
  }

  write.csv(stat_spector, file = paste(out, "SUMMARY_STAT_metric.csv", sep = ""),
            row.names = FALSE)
}

.save_merged <- function(res_v, out) {
    out_path <- paste(out, "results_bam_out_merged.csv", sep = "")
    write.csv(res_v, file = out_path, row.names = FALSE)
}

.read_id_assign <- function(id_path, f_head = FALSE) {
  list.name <- c("fs_bam", 'id_bam', 'gr_bam', 'baseline', 's_prep')
  for (i in 1:length(list.name)) {
    assign(list.name[i], NULL)
  }
  id_df <- dplyr::tbl_df(read.csv(id_path, sep = ",", header = f_head,
    stringsAsFactors = FALSE, strip.white = TRUE))
  colnames(id_df) <- list.name[1:ncol(id_df)]
  return(id_df)
}

dir.exists <- function(d) {
    de <- file.info(d)$isdir
    ifelse(is.na(de), FALSE, de)
}

.out_trailing_fix <- function(out_F) {
  tmp_out <- unlist(strsplit(out_F, split = ""))
  if (!(tmp_out[length(tmp_out)] == "/")) {
    out_F <- paste(out_F, "/", sep = "")
  }
  return(out_F)

}
