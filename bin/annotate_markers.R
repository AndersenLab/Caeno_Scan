#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
# Load Bim data , Note change based on where the files are
elegans_bim <- read.csv('test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP_ID", "CM", "BP", "A1", "A2"))
briggsae_bim <-  read.csv('test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP_ID", "CM", "BP", "A1", "A2"))
tropicalis_bim <- read.csv('test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP_ID", "CM", "BP", "A1", "A2"))
# Load GFF data, Note change based on where the files are
elegans_gff <- read.csv('test_data/c_elegans/genomes/PRJNA13758/WS283/csq/PRJNA13758.WS283.csq.chrI.gff3', sep='\t')
elegans_gff_filtered <- filter(elegans_gff, type == 'mRNA') # filtered for just mRNA
briggsae_gff <- read.csv('test_data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq/QX1410_nanopore.Feb2020.csq.chrI.gff3', sep='\t')
tropicalis_gff <- read.csv('test_data/c_tropicalis/genomes/NIC58_nanopore/June2021/csq/NIC58_nanopore.June2021.csq.chrI.gff3', sep='\t')
# Creating function
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
bim$Gene_ID <- NA
# Update columns based on overlaps
  if (length(overlaps) > 0) {
    bim$Intragenic[subjectHits(overlaps)] <- TRUE
    bim$Gene_ID[subjectHits(overlaps)] <- mRNA_gr$gene_id[queryHits(overlaps)]
  }
return(bim)
}
# Test function with elegans
annotated_elegans <- annotateSNPs(elegans_bim, elegans_gff_filtered, 100)
# Test function with briggsae
annotated_briggsae <- annotateSNPs(briggsae_bim, briggsae_gff, 100)
# Test function with tropicalis
annotated_tropicalis <- annotateSNPs(tropicalis_bim, tropicalis_gff, 100)
