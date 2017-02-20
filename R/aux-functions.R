#' @importFrom stringr str_detect ignore.case
#'
checkGenome <- function(genome) {
  if (genome %in% c("GRCh37", "hg19", "GRCh38", "hg38")) {
    if (str_detect(genome, ignore.case("GRCh"))) {
      if (str_detect(genome, "38")) {
        genome <- "hg38"
      } else if (str_detect(genome, "37")) {
        genome <- "hg19"
      }
    }

    message("Genome version selected: ", genome)
  } else {
    stop("This package only supports `GRCh37` and `GRCh38` genomes.\n",
       "  For any other genome version please provide a relevant bed file ",
       "using the `f_bed =` option")
  }

}
