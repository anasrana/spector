library(tidyverse)
library(stringr)

# ==========================================================================
# GIAB Data
# ==========================================================================

#
# Data source:
# https://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/analysis/GIAB_integration
#   GIAB v2.19_2
# --------------------------------------------------------------------------

url_giab <-
  str_c("https://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/analysis/",
    "GIAB_integration/union13callableMQonlymerged_addcert_nouncert_",
    "excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_",
    "noCNVs_v2.19_2mindatasets_5minYesNoRatio_AddRTGPlatGenConf_",
    "filtNISTclustergt9_RemNISTfilt_RemPartComp_RemRep_RemPartComp_v0.2.bed.gz")

# Download file and save in data-raw/
download.file(url_giab, "data-raw/giab.bed")

giab <-
  read.table("data-raw/giab.bed", col.names = c("chrom", "start", "end"),
    stringsAsFactors = FALSE) %>%
  tbl_df()

giab_10k <- giab %>%
  dplyr::filter(
    (end - start) > 10000, # filter out regions that are too large
    chrom %in% 1:22) %>%  # remove X-Chromosome
  dplyr::distinct(chrom, start, end)

# ==========================================================================
# GENOME DATA
# ==========================================================================

genome_version <- c("hg19", "hg38")


genome_size <-
lapply(genome_version, function(g_ver) {
  hg_db <- src_mysql(dbname = g_ver, host = "genome-mysql.cse.ucsc.edu",
                     user = "genome")


  hg_db %>%
    tbl("chromInfo") %>%
    select(chrom, size) %>%
    mutate(genome = g_ver) %>%
    collect() %>%
    dplyr::filter(str_detect(chrom, "chr[0-9]{1,2}\\b"))
  }) %>%
bind_rows()


genome_gap <-
lapply(genome_version, function(g_ver) {
  hg_db <- src_mysql(dbname = g_ver, host = "genome-mysql.cse.ucsc.edu",
                     user = "genome")
  hg_db %>%
    tbl("gap") %>%
    collect() %>%

    dplyr::filter(type  %in% c("telomere", "centromere")) %>%
    select(chrom, chromStart, chromEnd, type) %>%
    mutate(genome = g_ver) %>%
    collect() %>%
    dplyr::filter(str_detect(chrom, "chr[0-9]{1,2}\\b"))
  }) %>%
bind_rows()


# ==========================================================================
# Saving Data
# ==========================================================================

# Output giab_10k and genome data in R/sysdata.rda
devtools::use_data(giab_10k, genome_size, genome_gap, compress = "xz",
                   internal = TRUE, overwrite = TRUE)
