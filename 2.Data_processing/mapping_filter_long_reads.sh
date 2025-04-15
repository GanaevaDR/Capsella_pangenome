#!/bin/bash

for i in ./*.fastq.gz; do
    file=$(basename "$i" _trimmed.fastq.gz)
    minimap2 -ax map-pb ./reference/Co_Cr_reference.fasta "${file}_trimmed.fastq.gz" > ./"${file}.sam"

done

for i in ./*.sam; do
    file=$(basename "$i" .sam)

    samtools view -bSq 60 "${file}.sam" > "${file}.bam"
    samtools stats "${file}.bam" | grep ^SN | cut -f 2- > "${file}_stats.txt"

    samtools sort -@ 25 "${file}.bam" -o "${file}_sorted.bam"
    samtools index -@ 25 "${file}_sorted.bam"

done