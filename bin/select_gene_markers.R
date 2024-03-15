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

proc_dir = glue::glue("{out_dir}/snp_data")

# Check if the directory exists, if not create it
# if (!dir.exists(figure_dir)) {
#   dir.create(figure_dir)
# }

if (!dir.exists(proc_dir)) {
  dir.create(proc_dir)
}

params <- list(
  elegans_annotated = "c_elegans/ce.fullpop/gene_markers.tsv",
  briggsae_annotated = "c_briggsae/cb.fullpop/gene_markers.tsv",
  tropicalis_annotated = "c_tropicalis/ct.fullpop/gene_markers.tsv"
)


elegans_annotated <- read.csv(params$elegans_annotated, sep='\t' )
briggsae_annotated <-  read.csv(params$briggsae_annotated, sep='\t')
tropicalis_annotated <- read.csv(params$tropicalis_annotated, sep='\t')



select_snp <- function(df){
  # Group data by attribute
  grouped_data <- df %>% 
    group_by(gene_id) 
  
  # Filter data to keep only one SNP per attribute with lowest MAF
  filtered_data <- grouped_data %>%
    slice_min(MAF) %>%
    ungroup() %>%
    distinct(gene_id, .keep_all = TRUE)
  
  return(filtered_data)
  
  }

# Test function
# with elegans
ce_unique_groups <- length(unique(elegans_annotated$attribute))
selected_elegans <- select_snp(elegans_annotated)
ce_snps <- select(selected_elegans, marker)
ce_snps_maf <- select(selected_elegans, marker, MAF)
ce_selected_groups <- nrow(selected_elegans)

# with briggsae
cb_unique_groups <- length(unique(briggsae_annotated$attribute))
selected_briggsae <- select_snp(briggsae_annotated)
cb_snps <- select(selected_briggsae, marker)
cb_snps_maf <- select(selected_briggsae, marker, MAF)
cb_selected_groups <- nrow(selected_briggsae)

#with tropicalis
ct_unique_groups <- length(unique(tropicalis_annotated$attribute))
selected_tropicalis <- select_snp(tropicalis_annotated)
ct_snps <- select(selected_tropicalis, marker)
ct_snps_maf <- select(selected_tropicalis, marker, MAF)
ct_selected_groups <- nrow(selected_tropicalis)



print(
  glue::glue("There were {ce_unique_groups} attributes in the elegans data set, and 
  {ce_selected_groups} were selected.")
)


print(
  glue::glue("There were {cb_unique_groups} attributes in the briggsae data set, and 
  {cb_selected_groups} were selected.")
)


print(
  glue::glue("There were {ct_unique_groups} attributes in the tropicalis data set, and 
  {ct_selected_groups} were selected.")
)


# Sending out snp list
data.table::fwrite(ce_snps, glue::glue("{proc_dir}/{date}.snplist_elegans.tsv"))
data.table::fwrite(cb_snps, glue::glue("{proc_dir}/{date}.snplist__briggsae.tsv"))
data.table::fwrite(ct_snps, glue::glue("{proc_dir}/{date}.snplist__tropicalis.tsv"))



# Combine data for all species
ce_snps_maf$species <- "Elegans"
cb_snps_maf$species <- "Briggsae"
ct_snps_maf$species <- "Tropicalis"
all_data <- bind_rows(ce_snps_maf, cb_snps_maf, ct_snps_maf)
# Create histogram of the selected snps MAF and save it
c_allMAF_hist <- ggplot(all_data, aes(x = MAF, fill = species)) +
  geom_histogram(alpha = 0.7, bins = 30) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        title = element_text(size = 18, face = 'bold')) +
  theme_bw() +
  labs(title = "Histogram of MAF by Species",
       x = "MAF",
       y = "Count") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title =element_text(size=18, face='bold')) +
  facet_grid(. ~ species) +
  scale_fill_manual(values = c("Elegans" = "#DB6333", "Briggsae" = "#53886C", "Tropicalis" = "#0719BC"))

ggsave(
  glue::glue("{out_dir}/{date}.all_species_MAF.png"),
  plot = c_allMAF_hist,
  device = "png",
  width = 15,
  height = 9,
  units = "in",
  dpi = 300
)

