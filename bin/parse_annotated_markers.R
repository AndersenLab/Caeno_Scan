#Script to parse the annotated markers from the bedtools intersect command 
library(tidyverse)
#library(optparse)

args <- commandArgs(trailingOnly = TRUE)


#### Read Command line inputs
# path to the raw unannotated markers
markers_file <- args[1]
#markers_file <- "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim"

# path to the annoated .bim.bed.annotated file
anno_markers_file <- args[2]
#anno_markers_file <- "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim.bed.annotated" 

# path to the .freq plink file generated for these data
af_markers_file <- args[3]
#af_markers_file <- "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.frq"

# species id
sp_id <- args[4]
#sp_id <- "c_tropicalis"

ortho_gene_db_file <- "input_data/all_species/orthogroups/20240206_Orthogroups/masterOrthoDB_wAlias.tsv"

# Path to output
out_dir <- args[5]
#out_dir <- "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop"

###### 1. Parse the annotated markers file #####

#read in the annotation file

anno_df <- data.table::fread(
    anno_markers_file,
    col.names = c(
        "CHROM",
        "Start",
        "End",
        "marker",
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    )
    )

# Get the number of markers from the file output by bedtools 
n_markers_post_anno <- nrow(
    anno_df %>% 
        select(marker) %>% 
        unique()
    )

# Get the number of markers that was input into bedtools for a sanity check 
pre_bedtools_df <- data.table::fread(markers_file)  %>% 
    filter(V1 %in% c("1","2","3","4","5","6"))


n_markers_pre_anno <- nrow(pre_bedtools_df)

#get the SNPs that are within a feature of the GFF
overlapping_markers <- anno_df %>%
    dplyr::filter(seqid != ".") 

n_markers_overlap <- nrow(
    overlapping_markers %>% 
    select(marker) %>% 
    unique()
    )

#Now lets filter to just the overlaps between markers and gene features
overlapping_gene_markers <- overlapping_markers %>%
    dplyr::filter(type == "gene")

#Now we want to extract the gene name from the attributes column
#create a function that can be applied to the attributes column
get_name <- function(x) {
  id <- sub(";.*", "", x)
  transcript_id <- sub(".*:", "", id)
  return(transcript_id)
}

overlapping_gene_markers_clean <- overlapping_gene_markers %>%
  dplyr::mutate(gene_id = get_name(attributes))  %>% 
  select(CHROM, Start, End, marker, gene_id)

print(head(overlapping_gene_markers_clean))

n_gene_markers <- nrow(
    overlapping_gene_markers_clean %>% 
    select(marker) %>% 
    unique()
    ) 

#### 2. Add the orthogroup data
OG <- readr::read_tsv(ortho_gene_db_file) 
#colnames(OG) <- c("Orthogroup", "Briggsae", "Tropicalis", "Elegans")

# Now I want to get the orthogroup data in long format
if (sp_id == "c_elegans") {
    sp_ogs <- OG %>%
        select(Orthogroup, WB_id) %>% 
        rename("gene_id" = "WB_id")
} else if (sp_id == "c_briggsae") {
    sp_ogs <- OG %>%
        select(Orthogroup, QX1410)%>% 
        rename("gene_id" = "QX1410")
} else if (sp_id == "c_tropicalis") {
    sp_ogs <- OG %>%
        select(Orthogroup, NIC58)%>% 
        rename("gene_id" = "NIC58")
} else {
    print(glue::glue("Provided species id: {sp_id} is not found in {ortho_gene_db_file}"))
}

sp_ogs_long <- sp_ogs  %>% 
    separate_rows(gene_id, sep = ", ")


overlapping_gene_markers_clean_OG <- overlapping_gene_markers_clean %>%
    left_join(sp_ogs_long, by = c("gene_id" = "gene_id"))

#get just the snps that overlap with genes that are in an orthogroup (NA)
overlapping_gene_markers_orthogroups <- overlapping_gene_markers_clean_OG %>% 
    dplyr::filter(!is.na(Orthogroup) )

n_gene_markers_w_orthogroup <- nrow(overlapping_gene_markers_orthogroups)

# ct_ogs %>% 
#     select(Orthogroup, NIC58)  %>% 
#     separate_rows(NIC58, sep = ", ")  %>% 
#     head(n = 20)

#Now I want to join the orthogroup data to the gene markers





#### 3. Adding AF to the annotated data
af_data <- data.table::fread(af_markers_file)

#check the number of markers included in the allele frequency data
n_markers_af <- nrow(af_data)

#join the allele frequency data to the overlapping_gene_markers_clean_OG df
af_gene_markers <- overlapping_gene_markers_clean_OG  %>% 
    left_join(af_data, by = c("marker" = "SNP"))

write_tsv(
    af_gene_markers,
    glue::glue("{out_dir}/gene_markers.tsv")
    )


#### 4. Selecting representivitive SNP for each gene
#Function to select the snp with the highest
select_snp <- function(df){
  # Group data by attribute
  grouped_data <- df %>% 
    group_by(gene_id) 
  
  # Filter data to keep only one SNP per attribute with highest MAF
  filtered_data <- grouped_data %>%
    slice_max(MAF, with_ties = FALSE) %>%
    ungroup() 
  
  return(filtered_data)
  
  }

selected_snps <- select_snp(af_gene_markers)

plink_snps_in <- selected_snps  %>% 
    select(marker) %>% 
    unique() %>%
    pull(marker)

writeLines(
    plink_snps_in,
    glue::glue("{out_dir}/selected_snps.txt"),
)

#### 5. Output Variant Count Sanity Checks
var_counts <- data.frame(
    "n_markers_pre_ano" = c(n_markers_pre_anno),
    #the number of unique marker ids 
    "n_markers_post_anno" = c(n_markers_post_anno),
    #the number of markers in the allele frequency file
    "n_markers_af" = c(n_markers_af),
    #the number of unqique markers that overlap at least one feature in the gff
    "n_markers_overlap" = c(n_markers_overlap),
    #number of unique markers that overlap with a gene feature
    "n_markers_gene" = c(n_gene_markers),
    #number of markers that overlap with genes that are in OG
    "n_markers_gene_w_og" = c(n_gene_markers_w_orthogroup)
) 

print(var_counts)
