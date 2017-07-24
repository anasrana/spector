#' Las limits
#'
#' @param dat_df `tbl_df` object with LAS scores as a column.
#' @param las_lim `tbl_df` object containing upper and lower limits on LAS
#'        based on las_limits object contained in package.
#'
#' @return output is `tbl_df` object adding a column for region.status
#' @export
#'
#' @importFrom dplyr mutate if_else
las_status <- function(dat_df, las_lim) {

  dat_df %>%
    mutate(
      region.status = lasLimBin(las, las_lim) %>%
        if_else("abberation", "good"),
      region.status = if_else(las == 0 | is.na(las), as.character(NA),
        region.status)
      )

}
