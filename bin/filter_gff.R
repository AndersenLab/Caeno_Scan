#Code to filter the gffs to include only gene features for all three species 
library(tidyverse)
library(data.table)

ct_gff_file <- "input_data/c_tropicalis/annotations/c_tropicalis.NIC58_nanopore.June2021.csq.gff3"

#read in the gff file
read_gff <- function(gff_file_path){
  gff_df <- data.table::fread(
    gff_file_path,
    #set the column names
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  )
  #Get the total number of rows
    cat("Number of rows in the gff file: ", nrow(gff_df), "\n")
  
    #filter out the rows that are not genes
  gff_df <- gff_df %>%
    dplyr::filter(type == "gene")
    
    #Get the total number of rows
    cat("Number of rows in the genegff file: ", nrow(gff_df), "\n")
  
    #get the total number of unique rows
    cat("Number of unique rows in the gff file: ", nrow(unique(gff_df)), "\n")
  
  return(gff_df)
}

ct_gff_df <- read_gff(ct_gff_file)
head(ct_gff_df)

#save the filtered gff file
write.table(ct_gff_df, file = "input_data/c_tropicalis/annotations/c_tropicalis.NIC58_nanopore.June2021.csq.gene.gff3", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
