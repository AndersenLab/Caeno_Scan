#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(GenomicRanges)


# Load Bim data , Note change based on where the files are
elegans_bim <- read.csv('test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))
briggsae_bim <-  read.csv('test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))
tropicalis_bim <- read.csv('test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.bim', sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))

# Load GFF data, Note change based on where the files are
elegans_gff <- read.csv('test_data/c_elegans/genomes/PRJNA13758/WS283/csq/PRJNA13758.WS283.csq.chrI.gff3', sep='\t')
elegans_gff_filtered <- filter(elegans_gff, type == 'mRNA') # filtered for just mRNA

briggsae_gff <- read.csv('test_data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq/QX1410_nanopore.Feb2020.csq.chrI.gff3', sep='\t')
briggsae_gff_filtered <- filter(briggsae_gff, type == 'mRNA') # filtered for just mRNA

tropicalis_gff <- read.csv('test_data/c_tropicalis/genomes/NIC58_nanopore/June2021/csq/NIC58_nanopore.June2021.csq.chrI.gff3', sep='\t')
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
OG <- readr::read_tsv('input_data/all_species/orthogroups/20240206_Orthogroups/masterOrthoDB.tsv') 
colnames(OG) <- c("Orthogroup", "Briggsae", "Tropicalis", "Elegans")
cpOG <- OG
OG = unite(OG, Gene_ID, c("Briggsae", "Tropicalis", "Elegans"), sep = " ") %>% 
  separate_rows(Gene_ID, sep = ",") %>% 
  separate_rows(Gene_ID, sep = " ")
OG$Gene_ID <- sub("Transcript_", "", OG$Gene_ID)
OG$Gene_ID <- sub("transcript_", "", OG$Gene_ID) 


# function to add OGs 
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

AF_ce <- read.table('test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.chr1.frq', header = TRUE)
AF_cb <- read.table('test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.chr1.frq', header = TRUE)
AF_ct <- read.table('test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.chr1.frq', header = TRUE)


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

briggsae_hist<- ggplot(genes_briggsae, aes(x = Count)) + geom_histogram(fill = "deeppink", color = "black", binwidth = 1) +
  labs(title = "Histogram of Birggsae Counts w 100 bp buffer",
       x = "Number SNPs per Gene",
       y = "Frequency")
#print(briggsae_hist)
tropicalis_hist<- ggplot(genes_tropicalis, aes(x = Count)) + geom_histogram(fill = "blueviolet", color = "black", binwidth = 1) +
  labs(title = "Histogram of Tropicalis Counts w 100 bp buffer",
       x = "Number SNPs per Gene",
       y = "Frequency")
#print(tropicalis_hist)


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
  


# Venn Diagram
ce_one_one_one_var_ogs <- merge(one_one_one_og_var, filtered_all_elegans, by = "Orthogroup") 
ce_one_one_one_var_ogs <- subset(ce_one_one_one_var_ogs, Ratio == "1:1:1")

cb_one_one_one_var_ogs <- merge(one_one_one_og_var, filtered_all_briggsae, by = "Orthogroup") 
cb_one_one_one_var_ogs <- subset(ce_one_one_one_var_ogs, Ratio == "1:1:1")

ct_one_one_one_var_ogs <- merge(one_one_one_og_var, filtered_all_tropicalis, by = "Orthogroup") 
ct_one_one_one_var_ogs <- subset(ce_one_one_one_var_ogs, Ratio == "1:1:1")

library(VennDiagram)
#venn.diagram()









































