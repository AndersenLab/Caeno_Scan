This is code to generate test data for the marker set annotation project

# Filter the the GFF file to a single chromosome 
Start an interactive session on quest
```{bash}
srun --account=b1042 \
--partition=genomicsguest \
--time=02:00:00  \
--mem=30G \
--pty bash -l

module load R/4.1.1
R
```

Load libs
```{r}
library(tidyverse)
```

## C. elegans
```{r}
ce_gff_file <- "/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/csq/c_elegans.PRJNA13758.WS283.csq.gff3"

gff <- data.table::fread(
    ce_gff_file,
    col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
)

gff <- gff %>% 
    filter(chrom == "I") 

out_dir <- "test_data/c_elegans/genomes/PRJNA13758/WS283/csq"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

write.table(gff, 
    file = glue::glue("{out_dir}/PRJNA13758.WS283.csq.chrI.gff3"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
    )
rm(list=ls())
```

## C. briggsae
```{r}
cb_gff_file <- "/projects/b1059/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq/c_briggsae.QX1410_nanopore.Feb2020.csq.gff3"

gff <- data.table::fread(
    cb_gff_file,
    col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
)

gff <- gff %>% 
    filter(chrom == "I")

out_dir <- "test_data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

write.table(gff, 
    file = glue::glue("{out_dir}/QX1410_nanopore.Feb2020.csq.chrI.gff3"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
    )
rm(list=ls())
```

## C. tropicalis
```{r}
ct_gff_file <- "/projects/b1059/data/c_tropicalis/genomes/NIC58_nanopore/June2021/csq/c_tropicalis.NIC58_nanopore.June2021.csq.gff3"

gff <- data.table::fread(
    ct_gff_file,
    col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
)

gff <- gff %>% 
    filter(chrom == "I")

out_dir <- "test_data/c_tropicalis/genomes/NIC58_nanopore/June2021/csq"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

write.table(gff, 
    file = glue::glue("{out_dir}/NIC58_nanopore.June2021.csq.chrI.gff3"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
    )

rm(list=ls())
```