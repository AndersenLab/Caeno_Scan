#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(optparse)

#get date and time format as a variable YYYYMMDD_HHMM
date <- format(Sys.time(), "%Y%m%d_%H%M")

figure_dir = glue::glue("analysis/{date}_OG_SNPs_test/figures")

# Set up command line arguments
option_list = list(
  make_option(c("-e", "--elegans_bim"),  type="character"),
  make_option(c("-b", "--briggsae_bim"),  type="character"),
  make_option(c("-t", "--tropicalis_bim"),  type="character"),
  make_option(c("-E", "--elegans_gff"),  type="character"),
  make_option(c("-B", "--briggsae_gff"), type="character"),
  make_option(c("-T", "--tropicalis_gff"), type="character"),
  make_option(c("-f", "--elegans_freq"), type="character"),
  make_option(c("-g", "--briggsae_freq"), type="character"),
  make_option(c("-h", "--tropicalis_freq"), type="character"),
  make_option(c("-o", "--orthogroups"), type="character"),
  make_option(c("-O", "--output"), type="character")
)

# Parse the arguments
opt_parser = OptionParser(option_list=option_list)
params = parse_args(opt_parser)


# Load Bim data , Note change based on where the files are
elegans_bim <- read.csv(params$elegans_bim, sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))
briggsae_bim <-  read.csv(params$briggsae_bim, sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))
tropicalis_bim <- read.csv(params$tropicalis_bim, sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))

# Load GFF data, Note change based on where the files are
elegans_gff <- read.csv(params$elegans_gff, sep='\t')
elegans_gff_filtered <- filter(elegans_gff, type == 'mRNA') # filtered for just mRNA

briggsae_gff <- read.csv(params$briggsae_gff, sep='\t')
briggsae_gff_filtered <- filter(briggsae_gff, type == 'mRNA') # filtered for just mRNA

tropicalis_gff <- read.csv(params$tropicalis_gff, sep='\t')
tropicalis_gff_filtered <- filter(tropicalis_gff, type == 'mRNA') # filtered for just mRNA

# Creating function to annotate the snps 
annotateSNPs <- function(bim, gff, buffer){
  
# Convert PLINK chromosome to roman numeral to match GFF file
bim$chrom <- as.character(as.roman(bim$chrom))

# Create GRanges objects for SNPs and mRNA
snp_gr <- GRanges(seqnames = bim$chrom, ranges = IRanges(start = bim$BP - buffer, end = bim$BP + buffer))
mRNA_gr <- GRanges(seqnames = unique(gff$chrom), ranges = IRanges(start = gff$start - buffer, end = gff$end + buffer),
                   gene_id = gsub(".*transcript_id \"([^\"]+)\".*", "\\1", gff$attribute))

# Find overlapping SNPs and mRNA features
overlaps <- findOverlaps(snp_gr, mRNA_gr)

# Initialize result columns
bim$Intragenic <- FALSE
bim$attribute <- NA

# Update columns based on overlaps
  if (length(overlaps) > 0) {
    bim$Intragenic[subjectHits(overlaps)] <- TRUE
    bim$attribute[subjectHits(overlaps)] <- mRNA_gr$gene_id[queryHits(overlaps)]
  }

bim$Gene_ID <- sub(".*:(\\w+\\.\\w+\\.\\w+);.*", "\\1", bim[,8])

return(bim)
}


# Test annotate function with elegans
annotated_elegans <- annotateSNPs(elegans_bim, elegans_gff_filtered, 0)

# Test function with briggsae
annotated_briggsae <- annotateSNPs(briggsae_bim, briggsae_gff_filtered, 0)

# Test function with tropicalis
annotated_tropicalis <- annotateSNPs(tropicalis_bim, tropicalis_gff_filtered, 0)

# Read in average allele frequency data

AF_ce <- read.table(params$elegans_freq, header = TRUE)
AF_cb <- read.table(params$briggsae_freq, header = TRUE)
AF_ct <- read.table(params$tropicalis_freq, header = TRUE)


# create function to add MAFs
add_MAFs <- function(abim, allele_df) {
  # Filter AF frame to join
  filt_AF <- select(allele_df, SNP, MAF)
  
  # Join Data
  merged_data <- merge(abim, filt_AF, by = "SNP")
  
  return(merged_data)
}

# test function with elegans
all_elegans <- add_MAFs(OG_elegans, AF_ce)
# test function with briggsae
all_briggsae <- add_MAFs(OG_briggsae, AF_cb)
# test function with tropicalis
all_tropicalis <- add_MAFs(OG_tropicalis, AF_ct)

# Filter out NA values so just the GeneIDs

filtered_all_elegans <- all_elegans %>% filter(!is.na(Gene_ID))
filtered_all_briggsae <- all_briggsae %>% filter(!is.na(Gene_ID))
filtered_all_tropicalis <- all_tropicalis %>% filter(!is.na(Gene_ID))

# Count the occurrences of each Gene_ID
count_genes <- function(df) {
  gene_counts <- df %>% 
    group_by(Gene_ID, Orthogroup) %>% 
    summarise(count = n()) %>%
    filter(count > 1)
  # Create a new data frame with Gene_ID and their counts
  new_df <- data.frame(Gene_ID = gene_counts$Gene_ID, Count = gene_counts$count, Orthogroup = gene_counts$Orthogroup)

  return(new_df)
  
}


# test function
genes_elegans <- count_genes(filtered_all_elegans)
genes_briggsae <- count_genes(filtered_all_briggsae)
genes_tropicalis <- count_genes(filtered_all_tropicalis)


#Plot the counts
elegans_hist<- ggplot(genes_elegans, aes(x = Count)) + geom_histogram(fill = 'cyan', color = "black", binwidth = 1) +
  labs(title = "Histogram of Elegans Counts w 100 bp buffer",
       x = "Number SNPs per Gene",
       y = "Frequency")
#print(elegans_hist)

ggsave(
  glue::glue("{figure_dir}/{date}.elegans_hist.png"),
  plot = elegans_hist,
  device = "png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

briggsae_hist<- ggplot(genes_briggsae, aes(x = Count)) + geom_histogram(fill = "deeppink", color = "black", binwidth = 1) +
  labs(title = "Histogram of Birggsae Counts w 100 bp buffer",
       x = "Number SNPs per Gene",
       y = "Frequency")
#print(briggsae_hist)

ggsave(
  glue::glue("{figure_dir}/{date}.briggsae_hist.png"),
  plot = briggsae_hist,
  device = "png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

tropicalis_hist<- ggplot(genes_tropicalis, aes(x = Count)) + geom_histogram(fill = "blueviolet", color = "black", binwidth = 1) +
  labs(title = "Histogram of Tropicalis Counts w 100 bp buffer",
       x = "Number SNPs per Gene",
       y = "Frequency")
#print(tropicalis_hist)

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


## plot the data
all_counts_plot <- ggplot(all_counts, aes(x = Count)) + geom_histogram(fill = "cyan", color = "black", binwidth = 1) +
  labs(title = "Histogram of Counts w 100 bp buffer",
       x = "Number SNPs per Gene",
       y = "Frequency") +
  facet_grid(Species ~ ., scales = "fixed")

#print(all_counts_plot)

ggsave(
  glue::glue("{figure_dir}/{date}.all_counts_plot.png"),
  plot = all_counts_plot,
  device = "png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)






