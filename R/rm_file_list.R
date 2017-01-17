spectorFile <- function(f_bam,
                          id_bam = NULL,
                          srm_bam = NULL,
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
  srm_df$grp <- srm_bam
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

spectorList <- function(fs_bam, id_v, grp_v, s_v, out_F,
                          n_core = 1, save_out, ...) {
# function to be run inside of mclappy
  f_idx <- 1:length(fs_bam)
  spector_l <- parallel::mclapply(X = f_idx,
    function(idx) {
                            message(paste("Running file:", fs_bam[idx]))
                            spectorFile(
                              f_bam = fs_bam[idx],
                              id_bam = id_v[idx],
                              srm_bam = grp_v[idx],
                              s_prep = NULL,
                              out_F = out_F,
                              save_out = save_out,
                              ...)},
                                  mc.cores = n_core)
  srm_df <- dplyr::bind_rows(spector_l)
  return(srm_df)
}

