spectorFile <- function(f_bam,
                          id_bam = NULL,
                          s_prep = NULL,
                          out_F = NULL,
                          save_out,
                          ...) {

  message(paste("File:", f_bam, "\n=>\n"))

  add.args <- list(...)

  srm_df <- spector_metric(f_bam = f_bam, ...)
  if (is.null(id_bam)) {
    id_bam <- gsub(".bam", "", x = basename(f_bam))
  }

  srm_df$id_bam <- id_bam
  srm_df$prep <- s_prep


  if (is.null(out_F) & save_out) {
    warning("Output folder not specified, output will not be saved",
      call. = FALSE, domain = 'spector_run')
  } else if (save_out) {
    out_path <- paste(out_F, id_bam, "_out.csv", sep = "")
    write.csv(srm_df, file = out_path, row.names = FALSE)
    message(paste("Completed, output written to:", out_path))
  }

  return(srm_df)
}

#' Run spector for a list of bam file.
#'
#' @param fs_bam
#' @param id_v
#' @param grp_v
#' @param s_v
#' @param out_F
#' @param file_cores
#' @param save_out
#'
#' @return
#'
#' @importFrom parallel mclapply
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_c
#'
spectorList <- function(fs_bam, id_v, s_v, out_F,
                          file_cores = 1, save_out, ...) {
  # function to be run inside of mclappy
  f_idx <- 1:length(fs_bam)

  srm_df <-
  mclapply(X = f_idx, function(idx) {
    message(str_c("Running file: ", fs_bam[idx]))
      spectorFile(f_bam = fs_bam[idx], id_bam = id_v[idx],
        s_prep = s_v[idx], out_F = out_F, save_out = save_out, ...)},
    mc.cores = file_cores) %>%
  bind_rows()

  return(srm_df)
}

