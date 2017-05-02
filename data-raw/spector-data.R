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

url_giab_hg19 <-
  str_c("ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/",
    "NA12878_HG001/NISTv3.3.2/GRCh37/",
    "HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID",
    "_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed")

url_giab_hg38 <-
  str_c("ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/",
    "NA12878_HG001/NISTv3.3.2/GRCh38/",
    "HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID",
    "_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed")



# Download file and save in data-raw/
download.file(url_giab_hg19, "data-raw/giab_hg19.bed")
download.file(url_giab_hg38, "data-raw/giab_hg38.bed")

giab_hg19 <-
  read.table("data-raw/giab_hg19.bed", col.names = c("chrom", "start", "end"),
    stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  mutate(genome = "hg19")

giab_hg38 <-
  read.table("data-raw/giab_hg38.bed", col.names = c("chrom", "start", "end"),
    stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  mutate(genome = "hg38")

giab <-
  read.table("data-raw/giab.bed", col.names = c("chrom", "start", "end"),
    stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  mutate(genome = "hg38")

giab_10k <-
  bind_rows(giab_hg19, giab_hg38) %>%
  # remove X-Chromosome
  dplyr::filter(chrom %in% c(1:22, str_c("chr", 1:22))) %>%
  # filter out regions that are too large
  dplyr::filter((end - start) > 2^13) %>%
  group_by(genome) %>%
  dplyr::distinct(chrom, start, end, .keep_all = TRUE)

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
# Baseline data
# ==========================================================================

#
# Pre computed baselines for 2^13 - 2^17 for comparison and computing GIM
# --------------------------------------------------------------------------

base_gim <- read_csv("data-raw/baseline.csv", col_type = "ddi")


# ==========================================================================
# Saving Data
# ==========================================================================

# Output giab_10k and genome data in R/sysdata.rda
devtools::use_data(giab_10k, genome_size, genome_gap, base_gim, compress = "xz",
                   internal = TRUE, overwrite = TRUE)
