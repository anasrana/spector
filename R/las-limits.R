#' Las limits
#'
#' @param dat_df
#' @param las_lim
#'
#' @return
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
