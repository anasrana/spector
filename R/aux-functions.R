#' @importFrom stringr str_detect
#'
checkGenome <- function(genome) {
#ToDo check for implementation to allow for different case and mixed case
  if (str_detect(genome, "38|19|37")) {
    if (str_detect(genome, "38")) {
      genome <- "hg38"
    } else if (str_detect(genome, "37|19")) {
      genome <- "hg19"
    }
  }

  if (genome %in% c("hg19", "hg38")) {
    message("Genome version selected: ", genome)
  } else {
    stop("This package only supports `GRCh37` and `GRCh38` genomes.\n",
       "  For any other genome version please provide a relevant bed file ",
       "using the `f_bed =` option")
  }

}
