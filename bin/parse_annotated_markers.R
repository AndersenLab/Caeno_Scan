#Script to parse the annotated markers from the bedtools intersect command 
library(tidyverse)

#### Read Command line inputs
markers_file <- "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim"

anno_markers_file <- "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim.bed.annotated" 

af_markers_file <- ""

sp_id <- "c_tropicalis"

ortho_gene_db_file <- "input_data/all_species/orthogroups/20240206_Orthogroups/masterOrthoDB_wAlias.tsv"

#### 1. Parse the annotated markers file

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

n_markers_overlap <- nrow(overlapping_markers)

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

n_gene_markers <- nrow(overlapping_gene_markers_clean) 

#### 2. Add the orthogroup data
OG <- readr::read_tsv(ortho_gene_db_file) 
#colnames(OG) <- c("Orthogroup", "Briggsae", "Tropicalis", "Elegans")

# Now I want to get the orthogroup data in long format
if (sp_id == "c_elegans") {
    sp_ogs <- OG %>%
        select(Orthogroup, N2) %>% 
        rename("gene_id" = "N2")
}
if (sp_id == "c_briggsae") {
    sp_ogs <- OG %>%
        select(Orthogroup, QX1410)%>% 
        rename("gene_id" = "QX1410")
}
if (sp_id == "c_tropicalis") {
    sp_ogs <- OG %>%
        select(Orthogroup, NIC58)%>% 
        rename("gene_id" = "NIC58")
} else {
    print(glue::glue("Provided species id: {sp_id} is not found in {ortho_gene_db_file}"))
}


overlapping_gene_markers_clean_OG <- overlapping_gene_markers_clean %>%
    left_join(ct_ogs, by = c("gene_id" = "gene_id"))

# ct_ogs %>% 
#     select(Orthogroup, NIC58)  %>% 
#     separate_rows(NIC58, sep = ", ")  %>% 
#     head(n = 20)

#Now I want to join the orthogroup data to the gene markers


write_tsv(
    overlapping_gene_markers_clean_OG,
    "proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim.bed.annotated.gene_markers.tsv"
    )


#### 3. Adding AF to the annotated data


#### 4. Output Variant Count Sanity Checks
var_counts <- data.frame(
    "n_markers_pre_ano" = c(n_markers_pre_anno),
    "n_markers_post_anno" = c(n_markers_post_anno)
) 
