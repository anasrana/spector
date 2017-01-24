#' Get path for samples provided in the spector package.
#'
#' @param file file name
#'
#' @return path to \code{file}
#' @export
#'
#' @examples
#' spector_sample("sample1.bam")
spector_sample <- function(file) {
  system.file("extdata", file, package = "spector", mustWork = TRUE)
}
