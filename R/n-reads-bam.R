.nReadsBam <- function(f.name, cmd) {
    if (file.exists(paste(f.name, ".bai", sep = ""))) {
      n.read <- dplyr::tbl_df(read.table(
                  text = system(paste(cmd, " idxstats '", f.name, "'", sep = ""),
                  intern = TRUE, wait = TRUE), as.is = TRUE))
    } else {
      warning(".bam file is not indexed. Indexing is time consuming.")
        system(paste("cd '", dirname(f.name), "'&& ", cmd, " index '",
          basename(f.name), "'", sep =""), wait = TRUE)
      message("\nIndexing completed")
      n.read <- dplyr::tbl_df(read.table(
                  text = system(paste(cmd, " idxstats '", f.name, "'", sep = ""),
                  intern = TRUE, wait = TRUE), as.is = TRUE))
    }

    colnames(n.read) <- c("ref.name", 'seq.length', 'n.mapped', 'n.umapped')
    chrom <- (n.read %>%
        dplyr::filter(n.mapped > 0) %>%
        dplyr::distinct(ref.name))$ref.name
    n.read <- n.read %>%
        dplyr::mutate(n.read = n.mapped + n.umapped) %>%
        dplyr::summarise(n.read = sum(n.read / 10^9, na.rm = TRUE) )

    return(list(n.read = n.read$n.read, chrom = chrom))
}
