#' nReadsBam
#'
#' @param bam_file
#'
#' @importFrom Rsamtools indexBam
#'
nReadsBam <- function(bam_file) {
    if (file.exists(paste0(bam_file, ".bai"))) {
      bam_stat <- readIdxstats(bam_file)
    } else {
      warning(
        ".bam file is not indexed. Indexing is time consuming for large files.")
      indexBam(bam_file)
      message("Indexing completed")
      bam_stat <- readIdxstats(bam_file)
    }
    return(bam_stat)
}

#' readIdxstats
#'
#' @param bam_file
#'
#'
#' @importFrom dplyr filter mutate select
#' @importFrom Rsamtools idxstatsBam
#' @importFrom tibble as_data_frame
#'
readIdxstats <- function(bam_file) {
  idxstatsBam(bam_file) %>%
        as_data_frame(stringsAsFactors = FALSE) %>%
        mutate(
          n_read = mapped + unmapped,
          chrom = as.character(seqnames)) %>%
        filter(n_read > 0) %>%
        select(chrom, n_read)
}
