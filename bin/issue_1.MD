
Error with the CB SNP annotation script `annotate_markers.R`

The error appears to be associated with this code:
    
    ```R
    bim$Intragenic[subjectHits(overlaps)] <- TRUE
    bim$Gene_ID[subjectHits(overlaps)] <- mRNA_gr$gene_id[queryHits(overlaps)]
    ```

The SNP file has 

The CB bim file has 64640 rows 

The CB bim file is converted into a `GRanges` object wwith 64640 ranges. No elements appear to be lost doing this conversion.

Next, we look for overlaps with the gff

    ```R
    # Find overlapping SNPs and mRNA features
    overlaps <- findOverlaps(snp_gr, mRNA_gr)
    ```

The find overlaps function looks for elements in the query that overlap with the subject. This produces a hits object with 26105 elements.

Hits object containing the index pairings for the overlapping elements.

    ```R
    >print(head(overlaps))
    Hits object with 6 hits and 0 metadata columns:
        queryHits subjectHits
        <integer>   <integer>
    [1]         1          13
    [2]         1        4875
    [3]         1        4880
    [4]         1        4881
    [5]         1        4882
    [6]         1        4883
    -------
    queryLength: 64640 / subjectLength: 78587
    ```

So it looks like the first SNP overlaps with multiple mRNA features (13 and 4875)

```R
print(snp_gr[1])
GRanges object with 1 range and 0 metadata columns:
      seqnames      ranges strand
         <Rle>   <IRanges>  <Rle>
  [1]        I 57536-57736      *
  -------
  seqinfo: 7 sequences from an unspecified genome; no seqlengths

print(mRNA_gr[13])
GRanges object with 1 range and 1 metadata column:
      seqnames      ranges strand |                gene_id
         <Rle>   <IRanges>  <Rle> |            <character>
  [1]        I 55788-59063      * | ID=transcript:QX1410..
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

print(mRNA_gr[4875])
GRanges object with 1 range and 1 metadata column:
      seqnames      ranges strand |                gene_id
         <Rle>   <IRanges>  <Rle> |            <character>
  [1]        I 55788-59063      * | ID=gene:QX1410.4953;..
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

  "ID=gene:QX1410.4953;Name=QX1410.4953;locus=gpdh-3;biotype=protein_coding;Alias=gpdh-3,QX1410.4953"


print(mRNA_gr[4880])
GRanges object with 1 range and 1 metadata column:
      seqnames      ranges strand |                gene_id
         <Rle>   <IRanges>  <Rle> |            <character>
  [1]        I 56199-57619      * | ID=QX1410.4953.1-exo..
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

    "ID=QX1410.4953.1-exon3;Parent=transcript:QX1410.4953.1"

```

The error may be cause by the fact that the mRNA features in the CB gff file can 

This produces an object with 26105 elements.

The `queryHits(overlaps)` is returning a vector that is 26105 long. 

I beleive it is being used to subset the `mRNA_gr$gene_id` vector which is 78587 long.
