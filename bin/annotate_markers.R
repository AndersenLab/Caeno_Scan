#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
# Load Bim data , Note change based on where the files are
elegans_bim <- read.csv('~/Desktop/Erik/Caeno_Scan/test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP_ID", "CM", "BP", "A1", "A2"))
briggsae_bim <-  read.csv('~/Desktop/Erik/Caeno_Scan/test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP_ID", "CM", "BP", "A1", "A2"))
tropicalis_bim <- read.csv('~/Desktop/Erik/Caeno_Scan/test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP_ID", "CM", "BP", "A1", "A2"))

# Load GFF data, Note change based on where the files are
elegans_gff <- read.csv('~/Desktop/Erik/Caeno_Scan/test_data/c_elegans/genomes/PRJNA13758/WS283/csq/PRJNA13758.WS283.csq.chrI.gff3', sep='\t')
elegans_gff_filtered <- filter(elegans_gff, type == 'mRNA') # filtered for just mRNA

briggsae_gff <- read.csv('~/Desktop/Erik/Caeno_Scan/test_data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq/QX1410_nanopore.Feb2020.csq.chrI.gff3', sep='\t')
briggsae_gff_filtered <- filter(briggsae_gff, type == 'mRNA') # filtered for just mRNA

tropicalis_gff <- read.csv('~/Desktop/Erik/Caeno_Scan/test_data/c_tropicalis/genomes/NIC58_nanopore/June2021/csq/NIC58_nanopore.June2021.csq.chrI.gff3', sep='\t')
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
annotated_elegans <- annotateSNPs(elegans_bim, elegans_gff_filtered, 100)

# Test function with briggsae
annotated_briggsae <- annotateSNPs(briggsae_bim, briggsae_gff_filtered, 100)

# Test function with tropicalis
annotated_tropicalis <- annotateSNPs(tropicalis_bim, tropicalis_gff_filtered, 100)

#Read in orthogroups
OG <- readr::read_tsv('~/Desktop/Erik/Caeno_Scan/input_data/all_species/orthogroups/20240206_Orthogroups/masterOrthoDB.tsv') 
colnames(OG) <- c("Orthogroup", "A", "B", "C")
OG = unite(OG, Gene_ID, c("A", "B","C"), sep = " ") %>% 
  separate_rows(Gene_ID, sep = ",") %>% 
  separate_rows(Gene_ID, sep = " ")
OG$Gene_ID <- sub("Transcript_", "", OG$Gene_ID)
OG$Gene_ID <- sub("transcript_", "", OG$Gene_ID) 


# function with OGs 
add_OG <- function(abim, ortho) {
  
# Join the two tables based on 'Gene_ID'
merged_data <- left_join(abim, ortho, by = "Gene_ID")

# Group by 'Gene_ID' and make 'OG_list' as a list of unique OG values
result <- merged_data %>%
  group_by(Gene_ID) %>%
  summarize(OG_list = list(unique(Orthogroup, na.rm = TRUE)))

# Merge the summarized result back to the original 
final_data <- left_join(abim, result, by = "Gene_ID")

return(merged_data)

}

# test function with elegans
OG_elegans <- add_OG(annotated_elegans, OG)
# test function with briggsae
OG_briggsae <- add_OG(annotated_briggsae, OG)
# test function with tropicalis
OG_tropicalis <- add_OG(annotated_tropicalis, OG)


# Read in average allele frequency data

#allele_F <- read.table('~/Desktop/Erik/Caeno_Scan/test_data/<sp_id>/<c*>.comp.map/*.frq)

