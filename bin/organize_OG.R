#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(optparse)
library(VennDiagram)

#get date and time format as a variable YYYYMMDD_HHMM
date <- format(Sys.time(), "%Y%m%d_%H%M")

# # Parse the arguments
# opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
# params = parse_args(opt_parser)

# # Set up test data
out_dir = glue::glue("analysis/{date}_test_Orthogroups")

params <- list( 

  elegans_annotated = "c_elegans/ce.fullpop/gene_markers.tsv",
  briggsae_annotated = "c_briggsae/cb.fullpop/gene_markers.tsv",
  tropicalis_annotated = "c_tropicalis/ct.fullpop/gene_markers.tsv",
  out_dir = out_dir 
  )

### Uses annotated bim files from annotate_markers.R script ######

# check if the directory exists, if not create it
if (!dir.exists(params$out_dir)) {
  dir.create(params$out_dir)
}
out_dir = params$out_dir

figure_dir = glue::glue("{out_dir}/figures")

# Check if the directory exists, if not create it
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
}


# #Read in orthogroups
# OG <- readr::read_tsv('input_data/all_species/orthogroups/20240206_Orthogroups/masterOrthoDB.tsv') 
# colnames(OG) <- c("Orthogroup", "Briggsae", "Tropicalis", "Elegans")
# cpOG <- OG
# OG <- mutate(OG, Gene_ID =paste(Briggsae, Tropicalis, Elegans)) %>% 
#   select(-Briggsae) %>%
#   select(-Tropicalis) %>%
#   select(-Elegans) %>% 
#   separate_rows(Gene_ID, sep = ",") %>% 
#   separate_rows(Gene_ID, sep = " ")
# OG$Gene_ID <- sub("Transcript_", "", OG$Gene_ID)
# OG$Gene_ID <- sub("transcript_", "", OG$Gene_ID) 
# 
# # function to add OGs 
# add_OG <- function(abim, ortho) {
#   # Join the two tables based on 'Gene_ID'
#   merged_data <- left_join(abim, ortho, by = "Gene_ID")
#   # Group by 'Gene_ID' and make 'OG_list' as a list of unique OG values
#   result <- merged_data %>%
#     group_by(Gene_ID) %>%
#     summarize(OG_list = list(unique(Orthogroup, na.rm = TRUE)))
#   # Merge the summarized result back to the original 
#   final_data <- left_join(abim, result, by = "Gene_ID")
#   return(merged_data)
# }
# test function with elegans
# OG_elegans <- add_OG(annotated_elegans, OG)
# # test function with briggsae
# OG_briggsae <- add_OG(annotated_briggsae, OG)
# # test function with tropicalis
# OG_tropicalis <- add_OG(annotated_tropicalis, OG)

# create 1:1:1 OGs
# Function to determine the label based on count
get_ratio <- function(x) {
  if (is.na(x)) {
    return("none")
  } else {
    count <- length(unlist(strsplit(x, ", ")))
    if (count == 1) {
      return("1")
    } else {
      return("many")
    }
  }
}

# Apply the function to each row and create new column
cpOG$Ratio <- paste0(
  sapply(cpOG$Briggsae, get_ratio), ":",
  sapply(cpOG$Tropicalis, get_ratio), ":",
  sapply(cpOG$Elegans, get_ratio)
)

OG_structure <- cpOG[, c("Orthogroup", "Ratio")]
one_one_one_og_var <- subset(OG_structure, Ratio == "1:1:1")

# Load the annotated tsv files
ce_annotated <- readr::read_tsv(params$elegans_annotated)
cb_annotated <- readr::read_tsv(params$briggsae_annotated)
ct_annotated <- readr::read_tsv(params$tropicalis_annotated)

# Venn Diagram
ce_one_one_one_var_ogs <- merge(one_one_one_og_var, ce_bim, by = "Orthogroup")  %>%
  pull(Orthogroup)

cb_one_one_one_var_ogs <- merge(one_one_one_og_var, cb_bim, by = "Orthogroup") %>%
  pull(Orthogroup)

ct_one_one_one_var_ogs <- merge(one_one_one_og_var, ct_bim, by = "Orthogroup") %>%
  pull(Orthogroup)

venn.diagram(
  x = list(ce_one_one_one_var_ogs, cb_one_one_one_var_ogs, ct_one_one_one_var_ogs),
  category.names = c("C. elegans" , "C. briggsae " , "C. tropicalis"),
  filename = glue::glue("{figure_dir}/{date}_overlap_1_1_1_var_ogs.png"),
  output=TRUE,
  imagetype="png"
)
