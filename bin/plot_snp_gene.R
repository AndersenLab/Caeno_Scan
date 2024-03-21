#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(optparse)


#get date and time format as a variable YYYYMMDD_HHMM
date <- format(Sys.time(), "%Y%m%d_%H%M")

# Set up test data
out_dir = glue::glue("analysis/{date}_test_marker_gene")

params <- list( 
  elegans_annotated = "c_elegans/ce.fullpop/gene_markers.tsv",
  briggsae_annotated = "c_briggsae/cb.fullpop/gene_markers.tsv",
  tropicalis_annotated = "c_tropicalis/ct.fullpop/gene_markers.tsv",
  out_dir = out_dir 
)

# check if the directory exists, if not create it
if (!dir.exists(params$out_dir)) {
  dir.create(params$out_dir)
}
out_dir = params$out_dir

## Create sub directories for figures and processed data
figure_dir = glue::glue("{out_dir}/figures")
proc_dir = glue::glue("{out_dir}/proc_data")

# Check if the directory exists, if not create it
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
}

if (!dir.exists(proc_dir)) {
  dir.create(proc_dir)
}

# Load the annotated tsv files
ce_annotated <- readr::read_tsv(params$elegans_annotated)
cb_annotated <- readr::read_tsv(params$briggsae_annotated)
ct_annotated <- readr::read_tsv(params$tropicalis_annotated)



# Count the occurrences of each Gene_ID
count_genes <- function(df) {
  gene_counts <- df %>% 
    group_by(gene_id, Orthogroup) %>% 
    summarise(count = n()) %>%
    filter(count > 1)
  # Create a new data frame with Gene_ID and their counts
  new_df <- data.frame(Gene_ID = gene_counts$gene_id, Count = gene_counts$count, Orthogroup = gene_counts$Orthogroup)
  
  return(new_df)
  
}

# test function
genes_elegans <- count_genes(ce_annotated)
genes_briggsae <- count_genes(cb_annotated)
genes_tropicalis <- count_genes(ct_annotated)
#Plot the counts
elegans_hist<- ggplot(genes_elegans, aes(x = Count)) + geom_histogram(fill = "#DB6333", color = "black", binwidth = 1) +
  labs(title = "Histogram of Elegans Counts",
       x = "Number SNPs per Gene",
       y = "Frequency")

ggsave(
  glue::glue("{figure_dir}/{date}.elegans_hist.png"),
  plot = elegans_hist,
  device = "png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

briggsae_hist<- ggplot(genes_briggsae, aes(x = Count)) + geom_histogram(fill = "#53886C", color = "black", binwidth = 10) +
  labs(title = "Histogram of Birggsae Counts",
       x = "Number SNPs per Gene",
       y = "Frequency")


ggsave(
  glue::glue("{figure_dir}/{date}.briggsae_hist.png"),
  plot = briggsae_hist,
  device = "png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

tropicalis_hist<- ggplot(genes_tropicalis, aes(x = Count)) + geom_histogram(fill = "#0719BC", color = "black", binwidth = 1) +
  labs(title = "Histogram of Tropicalis Counts",
       x = "Number SNPs per Gene",
       y = "Frequency")


ggsave(
  glue::glue("{figure_dir}/{date}.tropicalis_hist.png"),
  plot = tropicalis_hist,
  device = "png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

# Plot on a shared axis using geom_facet_grid

## add a species column to each data frame
genes_elegans$Species <- "Elegans"
genes_briggsae$Species <- "Briggsae"
genes_tropicalis$Species <- "Tropicalis"

## combine all the data
all_counts <- rbind(genes_elegans, genes_briggsae, genes_tropicalis)

library(ggbreak) 
## plot the data
all_counts_plot <- ggplot(all_counts, aes(x = Count)) + geom_histogram(binwidth = 1) +
  labs(title = "Histogram of SNPs per Gene",
       x = "Number SNPs per Gene",
       y = "Frequency") +
  theme_bw() +
  facet_grid(Species ~ ., scales = "fixed") +
  scale_fill_manual(values = c("Elegans" = "#DB6333", "Briggsae" = "#53886C", "Tropicalis" = "#0719BC")) +
  ylim(0,500)+ #looking at lower frequencies, can modify this for what is needed 
  xlim(0,500)

all_counts_plot
ggsave(
  glue::glue("{figure_dir}/{date}.all_counts_plot.png"),
  plot = all_counts_plot,
  device = "png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)






