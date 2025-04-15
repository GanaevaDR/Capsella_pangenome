#!/bin/bash

for i in ./*_bwa_dna_pe_sorted.bam; do

    file=$(basename "$i" _bwa_dna_pe_sorted.bam)

    mosdepth -t 20 -b ../reference/Cori_chr_genes.bed -n "${file}_mean" ./"${file}_bwa_dna_pe_sorted.bam"

done