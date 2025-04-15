#!/bin/bash

for i in ./*_bwa_dna_pe_sorted.bam; do

    file=$(basename "$i" _bwa_dna_pe_sorted.bam)

    bedtools coverage -a ../reference/Cori_chr_genes.bed -b  ./"${file}_bwa_dna_pe_sorted.bam" > "res_${file}.txt"

done