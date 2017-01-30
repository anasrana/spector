#' spector: SEquence COverage Roughness
#'
#' spector provides an efficient and fast quality control for genomic data.
#'
#' It has three main goals:
#'
#' \itemize{
#' \item Identify regions in genomic data that are badly sequenced.
#' \item Determine sample quality.
#' \item Allow for easy comparison of quality between samples.
#' }
#'
#'
#' To learn more about spector, start with the vignettes:
#' `browseVignettes(package = "spector")`
#'
#' @seealso [spector_qc()] This is the main interface to run spector quality
#'          control on a bam file.
#'
#' @docType package
#' @name spector
#' @importFrom stats IQR cov end median sd setNames start
NULL
