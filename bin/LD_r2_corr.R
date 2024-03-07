#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(optparse)


#get date and time format as a variable YYYYMMDD_HHMM
date <- format(Sys.time(), "%Y%m%d_%H%M")

# Set up test data
out_dir = glue::glue("analysis/{date}_test_marker_gene")

params <- list( 
  
  ce_LD_r2 = "projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/ce.fullpop.ld.bin",
  cb_LD_r2 = "projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/ce.fullpop.ld.bin",
  ct_LD_r2 = "projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/ce.fullpop.ld.bin",
  out_dir = out_dir 
)

# check if the directory exists, if not create it
if (!dir.exists(params$out_dir)) {
  dir.create(params$out_dir)
}
out_dir = params$out_dir
