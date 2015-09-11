.spector_file <- function(f.bam,
                          id_bam = NULL,
                          srm_bam = NULL,
                          s_prep = NULL,
                          out_F = NULL,
                          stl_cmd = NULL,
                          ...) {

  message(paste("File:", f.bam, "\n=>\n"))

  add.args <- list(...)

  srm.df <- spector_metric(f.bam = f.bam, stl_cmd = stl_cmd, ...)
  if (is.null(id_bam)) {
    id_bam <- gsub(".bam", "", x = basename(f.bam))
  }

  srm.df$id <- rep(id_bam, nrow(srm.df))
  srm.df$grp <- rep(srm_bam, nrow(srm.df))
  srm.df$prep <- rep(s_prep, nrow(srm.df))


  if (is.null(out_F)) {
    warning("Output folder not specified, output will not be saved", call. = FALSE,
      domain = 'spector_run')
  } else {
    out_path <- paste(out_F, id_bam, "_out.csv", sep = "")
    write.csv(srm.df, file = out_path, row.names = FALSE)
    message(paste("Completed, output written to:", out_path))
  }

  return(srm.df)
}

.spector_list <- function(fs_bam, id_v, grp_v, s_v, out_F,
                          n_core = 1, stl_cmd, ...) {
# function to be run inside of mclappy
  f_idx <- 1:length(fs_bam)
  spector_l <- parallel::mclapply(X = f_idx, function(idx) {
                            message(paste("Running file:", fs_bam[idx]))
                            .spector_file(f.bam = fs_bam[idx], id_bam = id_v[idx],
                                  srm_bam = grp_v[idx], s_prep = NULL, out_F = out_F,
                                  stl_cmd = stl_cmd, ...)},
                                  mc.cores = n_core)
  srm.df <- dplyr::bind_rows(spector_l)
  return(srm.df)
}

