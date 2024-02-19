#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(optparse)
library(data.table)

#chrom = 'I'


#get date and time format as a variable YYYYMMDD_HHMM
date <- format(Sys.time(), "%Y%m%d_%H%M")

# Create the data folder in the analysis folder
out_dir = glue::glue("analysis/{date}_{chrom}_OG_SNPs_test")

# check if the directory exists, if not create it
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

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

# # Set up command line arguments
# option_list = list(
#   make_option(c("-e", "--elegans_bim"),  type="character"),
#   make_option(c("-b", "--briggsae_bim"),  type="character"),
#   make_option(c("-t", "--tropicalis_bim"),  type="character"),
#   make_option(c("-E", "--elegans_gff"),  type="character"),
#   make_option(c("-B", "--briggsae_gff"), type="character"),
#   make_option(c("-T", "--tropicalis_gff"), type="character"),
#   make_option(c("-f", "--elegans_freq"), type="character"),
#   make_option(c("-g", "--briggsae_freq"), type="character"),
#   make_option(c("-h", "--tropicalis_freq"), type="character")
# )

# # Parse the arguments
# opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
# params = parse_args(opt_parser)


# Set up inputs for troubleshooting error with large data set
# Set up inputs for troubleshooting error with large data set
params <- list(
  elegans_bim = "test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.bim",
  briggsae_bim = "test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.bim",
  tropicalis_bim = "test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.bim",

  elegans_gff = "test_data/c_elegans/genomes/PRJNA13758/WS283/csq/PRJNA13758.WS283.csq.chrI.gff3",
  briggsae_gff = "test_data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq/QX1410_nanopore.Feb2020.csq.chrI.gff3",
  tropicalis_gff = "test_data/c_tropicalis/genomes/NIC58_nanopore/June2021/csq/NIC58_nanopore.June2021.csq.chrI.gff3",

  elegans_freq = "test_data/c_elegans/ce.comp.map/ce.comp.map_0.05.chr1.frq",
  briggsae_freq = "test_data/c_briggsae/cb.comp.map/cb.comp.map_0.05.chr1.frq",
  tropicalis_freq = "test_data/c_tropicalis/ct.comp.map/ct.comp.map_0.05.chr1.frq"
)

# # Set up inputs for troubleshooting error with large data set
# params <- list(
#   elegans_bim = "/projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/20240212_fullpopulation_simfiles_noLD_0.00/c_elegans/ce_fullpop/ce_fullpop_0.00.bim",
#   briggsae_bim = "/projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/20240212_fullpopulation_simfiles_noLD_0.00/c_briggsae/cb_fullpop/cb_fullpop_0.00.bim",
#   tropicalis_bim = "/projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/20240212_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct_fullpop/ct_fullpop_0.00.bim",
#   
#   elegans_gff = "/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/csq/c_elegans.PRJNA13758.WS283.csq.gff3",
#   briggsae_gff = "/projects/b1059/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq/c_briggsae.QX1410_nanopore.Feb2020.csq.gff3",
#   tropicalis_gff = "/projects/b1059/data/c_tropicalis/genomes/NIC58_nanopore/June2021/csq/c_tropicalis.NIC58_nanopore.June2021.csq.gff3",
#   
#   elegans_freq = "/projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/20240212_fullpopulation_simfiles_noLD_0.00/c_elegans/ce_fullpop/ce_fullpop_0.00.frq",
#   briggsae_freq = "/projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/20240212_fullpopulation_simfiles_noLD_0.00/c_briggsae/cb_fullpop/cb_fullpop_0.00.frq",
#   tropicalis_freq= "/projects/b1059/projects/Ryan/ortholog_sims/Caeno_Scan/20240212_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct_fullpop/ct_fullpop_0.00.frq"
# )


# Load Bim data , Note change based on where the files are
elegans_bim <- read.csv(params$elegans_bim, sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))
briggsae_bim <-  read.csv(params$briggsae_bim, sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))
tropicalis_bim <- read.csv(params$tropicalis_bim, sep='\t', header = FALSE, col.names = c("chrom", "SNP", "CM", "BP", "A1", "A2"))

# Load GFF data, Note change based on where the files are

## Testing transcript ID parsing from attributes column
# test_attributes_string <- "ID=transcript:Y74C9A.3.1;Parent=gene:WBGene00022277;Name=Y74C9A.3.1;wormpep=CE28146;locus=homt-1;uniprot_id=Q9N4D9;biotype=protein_coding"

# #extract everything before the first semicolon
# id <- sub(";.*", "", test_attributes_string)

# #pull out the transcript id from the id string
# transcript_id <- sub(".*:", "", id)

#create a function that can be applied to the attributes column
get_name <- function(x) {
  id <- sub(";.*", "", x)
  transcript_id <- sub(".*:", "", id)
  return(transcript_id)
}

elegans_gff_filtered <- data.table::fread(params$elegans_gff,
                                          sep='\t',
                                          header = FALSE, 
                                          col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "transcript_id"),
                                          # colClasses = c("character", "character", "character", "integer", "integer", "character", "character", "character", "character")
) %>% 
  dplyr::filter(type == 'mRNA') %>%
  dplyr::mutate(transcript_id = get_name(transcript_id)) %>% 
  #convert stop and end to integer
  mutate(start = as.integer(start), end = as.integer(end))

# Check if any of the transcript ids are NA
elegan_missing_id <- nrow(elegans_gff_filtered %>% filter(is.na(transcript_id)))


briggsae_gff_filtered <- data.table::fread(params$briggsae_gff,
                                           sep='\t',
                                           header = FALSE, 
                                           col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "transcript_id"),
                                           # colClasses = c("character", "character", "character", "integer", "integer", "character", "character", "character", "character")
) %>% 
  dplyr::filter(type == 'mRNA') %>% 
  dplyr::mutate(transcript_id = get_name(transcript_id)) %>% 
  #convert stop and end to integer
  mutate(start = as.integer(start), end = as.integer(end))

briggsae_missing_id <- nrow(briggsae_gff_filtered %>% filter(is.na(transcript_id)))

tropicalis_gff_filtered <- data.table::fread(params$tropicalis_gff,
                                             sep='\t',
                                             header = FALSE, 
                                             col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "transcript_id"),
                                             # colClasses = c("character", "character", "character", "integer", "integer", "character", "character", "character", "character")
) %>% 
  dplyr::filter(type == 'mRNA') %>%
  dplyr::mutate(transcript_id = get_name(transcript_id)) %>% 
  #convert stop and end to integer
  mutate(start = as.integer(start), end = as.integer(end))

tropicalis_missing_id <- nrow(tropicalis_gff_filtered %>% filter(is.na(transcript_id)))

annotateSNPs_bedtools <- function(bim, gff, buffer) {
  # Convert PLINK chromosome to roman numeral to match GFF file
  bim$chrom <- as.character(as.roman(bim$chrom))
  
  # Create BED data frames for SNPs and mRNA
  snp_bed <- data.frame(chrom = bim$chrom, start = bim$BP - buffer, end = bim$BP + buffer)
  mRNA_bed <- data.frame(chrom = gff$chrom, start = gff$start - buffer, end = gff$end + buffer, transcript_id = gff$transcript_id)
  
  # Ensure both data frames have the same number of columns
  max_cols <- max(ncol(snp_bed), ncol(mRNA_bed))
  snp_bed <- cbind(snp_bed, matrix(NA, nrow = nrow(snp_bed), ncol = max_cols - ncol(snp_bed)))
  mRNA_bed <- cbind(mRNA_bed, matrix(NA, nrow = nrow(mRNA_bed), ncol = max_cols - ncol(mRNA_bed)))

  # Run bedtools intersect 
  overlaps <- system2(
    command = "bedtools",
    args = c("intersect", "-wa", "-wb", "-loj"),
    stdin = capture.output(write.table(rbind(snp_bed, mRNA_bed), sep = "\t", col.names = FALSE, row.names = FALSE)),
    stdout = TRUE,
    stderr = TRUE
  )
  
  # Read the output
  overlaps <- read.table(text = overlaps, header = FALSE, col.names = c("snp_chr", "snp_start", "snp_end", "snp_index", "mRNA_chr", "mRNA_start", "mRNA_end", "transcript_id"), stringsAsFactors = FALSE)
  
  # Initialize result columns
  bim$Intragenic <- FALSE
  bim$attribute <- NA
  
  # Update columns based on overlaps
  if (nrow(overlaps) > 0) {
    for (i in 1:nrow(overlaps)) {
      snp_index <- overlaps$snp_index[i]
      mRNA_id <- overlaps$transcript_id[i]
      bim$Intragenic[snp_index] <- TRUE
      bim$attribute[snp_index] <- mRNA_id
    }
  }
  
  return(bim)
}



# Test annotate function with elegans
annotated_elegans <- annotateSNPs(elegans_bim, elegans_gff_filtered, 1)

#check the output
n_snps_elegans <- nrow(elegans_bim)
n_snps_annotated_elegans <- nrow(annotated_elegans)
n_snps_intragenic_elegans <- nrow(annotated_elegans %>% filter(Intragenic == TRUE))
n_snps_attribute_id_elegans <- nrow(annotated_elegans %>% filter(!is.na(attribute)))
#n_snps_geneid_elegans <- nrow(annotated_elegans %>% filter(!is.na(Gene_ID)))

print(
  glue::glue("There were {n_snps_elegans} SNPs in the elegans data set, 
  {n_snps_annotated_elegans} were annotated, 
  {n_snps_intragenic_elegans} were intragenic,
  and {n_snps_attribute_id_elegans} had an attribute id associated with them"
             # and {n_snps_geneid_elegans} had a gene id associated with them"
             
  )
)

# check output for intragenic SNPs without gene id
intragenic_no_gene_id <- annotated_elegans %>% filter(Intragenic == TRUE, is.na(attribute))
print(
  glue::glue("There were {nrow(intragenic_no_gene_id)} intragenic SNPs without a gene id")
)


# Test function with briggsae
annotated_briggsae <- annotateSNPs(briggsae_bim, briggsae_gff_filtered, 0)

#check the output
n_snps_briggsae <- nrow(briggsae_bim)
n_snps_annotated_briggsae <- nrow(annotated_briggsae)
n_snps_intragenic_briggsae <- nrow(annotated_briggsae %>% filter(Intragenic == TRUE))
print(
  glue::glue("There were {n_snps_briggsae} SNPs in the briggsae data set, {n_snps_annotated_briggsae} were annotated, and {n_snps_intragenic_briggsae} were intragenic")
)

# Test function with tropicalis
annotated_tropicalis <- annotateSNPs(tropicalis_bim, tropicalis_gff_filtered,0)

#check the output
n_snps_tropicalis <- nrow(tropicalis_bim)
n_snps_annotated_tropicalis <- nrow(annotated_tropicalis)
n_snps_intragenic_tropicalis <- nrow(annotated_tropicalis %>% filter(Intragenic == TRUE))
print(
  glue::glue("There were {n_snps_tropicalis} SNPs in the tropicalis data set, {n_snps_annotated_tropicalis} were annotated, and {n_snps_intragenic_tropicalis} were intragenic")
)

# Read in average allele frequency data

AF_ce <- read.table(params$elegans_freq, header = TRUE)
AF_cb <- read.table(params$briggsae_freq, header = TRUE)
AF_ct <- read.table(params$tropicalis_freq, header = TRUE)

# check nrow read in for each species
n_af_ce <- nrow(AF_ce)
n_af_cb <- nrow(AF_cb)
n_af_ct <- nrow(AF_ct)

print(
  glue::glue("There were {n_af_ce} SNPs in the elegans allele frequency data set, {n_af_cb} in the briggsae allele frequency data set, and {n_af_ct} in the tropicalis allele frequency data set")
)

# create function to add MAFs
add_MAFs <- function(abim, allele_df) {
  # Filter AF frame to join
  filt_AF <- select(allele_df, SNP, MAF)
  
  # Join Data
  merged_data <- merge(abim, filt_AF, by = "SNP")
  
  return(merged_data)
}

#Read in orthogroups
OG <- readr::read_tsv('input_data/all_species/orthogroups/20240206_Orthogroups/masterOrthoDB.tsv') 
colnames(OG) <- c("Orthogroup", "Briggsae", "Tropicalis", "Elegans")
cpOG <- OG
OG = unite(OG, attribute, c("Briggsae", "Tropicalis", "Elegans"), sep = " ") %>% 
  separate_rows(attribute, sep = ",") %>% 
  separate_rows(attribute, sep = " ")
OG$attribute <- sub("Transcript_", "", OG$attribute)
OG$attribute <- sub("transcript_", "", OG$attribute) 


# function to add OGs 
add_OG <- function(abim, ortho) {
  
  # Join the two tables based on 'Gene_ID'
  merged_data <- left_join(abim, ortho, by = "attribute")
  
  # Group by 'Gene_ID' and make 'OG_list' as a list of unique OG values
  result <- merged_data %>%
    group_by(attribute) %>%
    summarize(OG_list = list(unique(Orthogroup, na.rm = TRUE)))
  
  # Merge the summarized result back to the original 
  final_data <- left_join(abim, result, by = "attribute")
  
  return(merged_data)
  
}

# test function with elegans
OG_elegans <- add_OG(annotated_elegans, OG)

# check the output
n_snps_elegans <- nrow(annotated_elegans)
n_snps_annotated_elegans <- nrow(OG_elegans)
n_snps_intragenic_elegans <- nrow(OG_elegans %>% filter(Intragenic == TRUE))

print(
  glue::glue("There were {n_snps_elegans} SNPs in the OG elegans dataframe, {n_snps_annotated_elegans} were annotated with the OG and stored in the OG elegans dataframe,
  and {n_snps_intragenic_elegans} were intragenic")
)

# test function with briggsae
OG_briggsae <- add_OG(annotated_briggsae, OG)

# check the output
n_snps_briggsae <- nrow(annotated_briggsae)
n_snps_annotated_briggsae <- nrow(OG_briggsae)
n_snps_intragenic_briggsae <- nrow(OG_briggsae %>% filter(Intragenic == TRUE))

print(
  glue::glue("There were {n_snps_briggsae} SNPs in the OG briggsae dataframe, {n_snps_annotated_briggsae} were annotated with the OG and stored in the OG briggsae dataframe,
  and {n_snps_intragenic_briggsae} were intragenic")
)

# test function with tropicalis
OG_tropicalis <- add_OG(annotated_tropicalis, OG)

# check the output
n_snps_tropicalis <- nrow(annotated_tropicalis)
n_snps_annotated_tropicalis <- nrow(OG_tropicalis)
n_snps_intragenic_tropicalis <- nrow(OG_tropicalis %>% filter(Intragenic == TRUE))

print(
  glue::glue("There were {n_snps_tropicalis} SNPs in the OG tropicalis dataframe, {n_snps_annotated_tropicalis} were annotated with the OG and stored in the OG tropicalis dataframe,
  and {n_snps_intragenic_tropicalis} were intragenic")
)
# test function with elegans
all_elegans <- add_MAFs(OG_elegans, AF_ce)

#check the output
n_snps_elegans <- nrow(OG_elegans)
n_snps_annotated_elegans <- nrow(all_elegans)
n_snps_intragenic_elegans <- nrow(all_elegans %>% filter(Intragenic == TRUE))
print(
  glue::glue("There were {n_snps_elegans} SNPs in the OG elegans dataframe, {n_snps_annotated_elegans} were annotated with the MAF and stored in the all elegans dataframe,
  and {n_snps_intragenic_elegans} were intragenic")
)

# test function with briggsae
all_briggsae <- add_MAFs(OG_briggsae, AF_cb)

#check the output
n_snps_briggsae <- nrow(OG_briggsae)
n_snps_annotated_briggsae <- nrow(all_briggsae)
n_snps_intragenic_briggsae <- nrow(all_briggsae %>% filter(Intragenic == TRUE))
print(
  glue::glue("There were {n_snps_briggsae} SNPs in the OG briggsae dataframe, {n_snps_annotated_briggsae} were annotated with the MAF and stored in the all briggsae dataframe,
  and {n_snps_intragenic_briggsae} were intragenic")
)

# test function with tropicalis
all_tropicalis <- add_MAFs(OG_tropicalis, AF_ct)

#check the output
n_snps_tropicalis <- nrow(OG_tropicalis)
n_snps_annotated_tropicalis <- nrow(all_tropicalis)
n_snps_intragenic_tropicalis <- nrow(all_tropicalis %>% filter(Intragenic == TRUE))

print(
  glue::glue("There were {n_snps_tropicalis} SNPs in the OG tropicalis dataframe, {n_snps_annotated_tropicalis} were annotated with the MAF and stored in the all tropicalis dataframe,
  and {n_snps_intragenic_tropicalis} were intragenic")
)

# Filter out NA values so just the GeneIDs

filtered_all_elegans <- all_elegans %>% filter(!is.na(attribute))
filtered_all_briggsae <- all_briggsae %>% filter(!is.na(attribute))
filtered_all_tropicalis <- all_tropicalis %>% filter(!is.na(attribute))

# Save the data for later use
data.table::fwrite(filtered_all_elegans, glue::glue("{proc_dir}/{date}_{chrom}.filtered_all_elegans.tsv"), sep = "\t")
data.table::fwrite(filtered_all_briggsae, glue::glue("{proc_dir}/{date}_{chrom}.filtered_all_briggsae.tsv"), sep = "\t")
data.table::fwrite(filtered_all_tropicalis, glue::glue("{proc_dir}/{date}_{chrom}.filtered_all_tropicalis.tsv"), sep = "\t")

# Count the occurrences of each Gene_ID
count_genes <- function(df) {
  gene_counts <- df %>% 
    group_by(attribute, Orthogroup) %>% 
    summarise(count = n()) %>%
    filter(count > 1)
  # Create a new data frame with Gene_ID and their counts
  new_df <- data.frame(attribute = gene_counts$attribute, Count = gene_counts$count, Orthogroup = gene_counts$Orthogroup)
  
  return(new_df)
  
}


# test function
genes_elegans <- count_genes(filtered_all_elegans)
genes_briggsae <- count_genes(filtered_all_briggsae)
genes_tropicalis <- count_genes(filtered_all_tropicalis)


#Plot the counts
elegans_hist<- ggplot(genes_elegans, aes(x = Count)) + geom_histogram(fill = 'cyan', color = "black", binwidth = 1) +
  labs(title = "Histogram of Elegans Counts",
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
  labs(title = "Histogram of Birggsae Counts",
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
  labs(title = "Histogram of Tropicalis Counts",
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



