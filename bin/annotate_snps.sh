#!/bin/bash

marker_bed=$1
gff=$2
out_file=$3
#marker_bed=proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim.bed
#ct_prot_gff=input_data/c_tropicalis/annotations/c_tropicalis.NIC58_nanopore.June2021.csq.gene.gff3
#full_gff=input_data/c_tropicalis/annotations/c_tropicalis.NIC58_nanopore.June2021.csq.gff3

#Get folder name containing the marker_bed


#bedtools intersect -a $ct_maker_bed -b $ct_prot_gff -wa -wb > proc_data/20240304_fullpopulation_simfiles_noLD_0.00/c_tropicalis/ct.fullpop/ct.fullpop_0.00.bim.bed.annotated

bedtools intersect -a $marker_bed -b $gff -wa -wb -loj > $out_file
