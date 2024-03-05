library(tidyverse)
#Script to convert the bim file to a format it can be used to check overlaps with bedtools
bim_file_path <- "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim"

bim_df <- data.table::fread(
    bim_file_path,
    col.names = c(
        "chrom",
        "SNP",
        "CM",
        "BP",
        "A1",
        "A2"
        )
    )

#function to convert the bim file to a bed file

bim_to_bed <- function(bim_df){
    bed_df <- bim_df %>%
        filter(chrom %in% c("1", "2", "3", "4", "5", "6")) %>%
        dplyr::mutate(
            start = BP,
            end = BP,
            name = SNP
            ) %>% 
        dplyr::select(chrom, start, end, name)

    #rename the chromosome to roman numerals
    chrom_key <- c(
        "1" = "I",
        "2" = "II",
        "3" = "III",
        "4" = "IV",
        "5" = "V",
        "6" = "X"
    )
    
    bed_df <- bed_df %>%
        dplyr::mutate(chrom = chrom_key[chrom])

    return(bed_df)
}

marker_bed <- bim_to_bed(bim_df) 

#save the bed file
## get the directory of the bim file
bim_file_path_dir <- dirname(bim_file_path)

## save the bed file in the same directory as the bim file
write.table(
    marker_bed, 
    file = paste0(bim_file_path_dir, "/", basename(bim_file_path), ".bed"), 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE,
    col.names = FALSE)
