library(tidyverse)
library(stringr)

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

# Output giab_10k data_frame in R/sysdata.rda
devtools::use_data(giab_10k, compress = "xz", internal = TRUE, overwrite = TRUE)
