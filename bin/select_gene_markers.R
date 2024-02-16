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

# Create the data folder in the analysis folder
out_dir = glue::glue("analysis/{date}_select_markers_test")

# check if the directory exists, if not create it
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

## Create sub directories for figures and processed data

proc_dir = glue::glue("{out_dir}/proc_data")

# Check if the directory exists, if not create it
# if (!dir.exists(figure_dir)) {
#   dir.create(figure_dir)
# }

if (!dir.exists(proc_dir)) {
  dir.create(proc_dir)
}

params <- list(
  elegans_annotated = "test_data/proc_data/20240216_0836.filtered_all_elegans.tsv",
  briggsae_annotated = "test_data/proc_data/20240216_0836.filtered_all_briggsae.tsv",
  tropicalis_annotated = "test_data/proc_data/20240216_0836.filtered_all_tropicalis.tsv"
)


elegans_annotated <- read.csv(params$elegans_annotated, sep='\t' )
briggsae_annotated <-  read.csv(params$briggsae_annotated, sep='\t')
tropicalis_annotated <- read.csv(params$tropicalis_annotated, sep='\t')



select_snp <- function(df){
  # Group data by attribute
  grouped_data <- df %>% 
    group_by(attribute)
  
  # Filter data to keep only one SNP per attribute with lowest MAF
  filtered_data <- grouped_data %>%
    slice_min(MAF) %>%
    ungroup() %>%
    distinct(attribute, .keep_all = TRUE)
  
  return(filtered_data)
  
  }


# Test function
selected_elegans <- select_snp(elegans_annotated)
selected_briggsae <- select_snp(briggsae_annotated)
selected_tropicalis <- select_snp(tropicalis_annotated)


