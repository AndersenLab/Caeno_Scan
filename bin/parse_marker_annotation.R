library(tidyverse)
library(data.table)

#parse the command line arguments
args <- commandArgs(trailingOnly = TRUE)
#sp_annotation_file <- args[1]
sp_annotation_file <- "/Volumes/Macintosh HD/Users/ryanmckeown/Desktop/Caeno_Scan/input_data/c_tropicalis/annotations/WI.20231201.strain-annotation.tsv"

#marker_file <- args[2]
marker_file <- "proc_data/ct_fullpop_0.00.bim"

#read in the annotation file
read_anno <- function(anno_file_path, chroms = c("I", "II", "III", "IV", "V", "X")) {
    chrom_key <- c(
        "I" = "1",
        "II" = "2",
        "III" = "3",
        "IV" = "4",
        "V" = "5",
        "X" = "6"
    )
    anno_df <- data.table::fread(anno_file_path,
        select = c(
            "CHROM",
            "POS",
            "CONSEQUENCE",
            "WORMBASE_ID",
            "GENE",
            "TRANSCRIPT",
            "REF",
            "ALT"
        )
    ) %>%
        dplyr::filter(CHROM %in% chroms) %>%
        # filter out any rows with an '@' symbol in the Consequece column
        dplyr::filter(!grepl("@", CONSEQUENCE)) %>%
        dplyr::mutate(CHROM = chrom_key[CHROM]) %>%
        dplyr::mutate(marker = paste0(CHROM, ":", POS)) %>%
        dplyr::select(-CHROM, -POS, -CONSEQUENCE)

    print(names(anno_df))
    # Get the total number of rows
    cat("Number of rows in the annotation file: ", nrow(anno_df), "\n")

    # get the total number of unique rows
    cat("Number of unique rows in the annotation file: ", nrow(unique(anno_df)), "\n")

    return(anno_df)
}


#read in the marker file
read_bim <- function(bim_file_path){
  bim_df <- data.table::fread(
    bim_file_path,
    #set the column names
    col.names = c("CHROM", "SNP", "CM", "BP", "A1", "A2")
    )

    #Get the total number of rows
    cat("Number of rows in the bim file: ", nrow(bim_df), "\n")

    #get the total number of unique rows
    cat("Number of unique rows in the bim file: ", nrow(unique(bim_df)), "\n")

  return(bim_df)

}


anno_df <- read_anno(sp_annotation_file) %>%
    dplyr::select(-TRANSCRIPT)%>%
    unique()

#check how many SNPs are in the annotation file with a unique WORMBASE_ID
cat("Number of markers in the annotation file with a unique WORMBASE_ID: ", sum(!duplicated(anno_df$marker)), "\n")

bim_df <- read_bim(marker_file)

merged_df <- dplyr::left_join(
    bim_df,
    unique(anno_df),
    by = c(
        "SNP" = "marker",
        "A1" = "REF",
        "A2" = "ALT"
        )
    )

#Get the total number of rows
cat("Number of rows in the merged file: ", nrow(merged_df), "\n")

#Get the number of markers that do not have an annotation
cat("Number of markers that do not have an annotation: ", sum(is.na(merged_df$WORMBASE_ID)), "\n")

#Get the number of markers that have a non-NA WORMBASE_ID
cat("Number of markers that have an annotation: ", sum(!is.na(merged_df$WORMBASE_ID)), "\n")
#Get the number of markers that have more than one annotation ** DOES NOT WORK **
#cat("Number of markers that have more than one annotation: ", sum(duplicated(merged_df$SNP)), "\n")

### There should not be any makres where NA is in genes
multi_anno_markers <- merged_df%>%
  dplyr::group_by(SNP)%>%
  dplyr::summarize(
    n = n(),
    #print the gene names
    genes = paste0(unique(WORMBASE_ID), collapse = ", ")
    )%>%
  dplyr::filter(n > 1)

# print a preview of the multi_anno_markers
print(head(multi_anno_markers))

cat("Number of markers that have more than one annotation: ", nrow(multi_anno_markers), "\n")

write.table(
    merged_df,
    file = "merged_anno.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
    )