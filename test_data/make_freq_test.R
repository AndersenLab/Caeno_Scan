# Purpose: Create a test frequency file for the C. elegans data
# These .frq files were created on QUEST but are too large to share on GitHub

library(tidyverse)

## C. elegans ## ----
ce_freqs <- data.table::fread("test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.frq") 

ce_chrom_freqs <- ce_freqs %>% 
    filter(CHR == "1")

data.table::fwrite(
    ce_chrom_freqs,
    "test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.chr1.frq",
    sep = "\t"
    )

## C. briggsae ## ----
cb_freqs <- data.table::fread("test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.frq")

cb_chrom_freqs <- cb_freqs %>% 
    filter(CHR == "1")

data.table::fwrite(
    cb_chrom_freqs,
    "test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.chr1.frq",
    sep = "\t"
    )

## C. tropicalis ## ----
ct_freqs <- data.table::fread("test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.frq")

ct_chrom_freqs <- ct_freqs %>% 
    filter(CHR == "1")
    
data.table::fwrite(
    ct_chrom_freqs,
    "test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.chr1.frq",
    sep = "\t"
    )
    