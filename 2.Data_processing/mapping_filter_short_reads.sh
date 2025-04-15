#!/bin/bash

for i in ./*1_trimmed.fastq.gz; do

    file=$(basename "$i" _1_trimmed.fastq.gz)

    bwa-mem2 mem -t 25 ./reference/bwa_index/Co_index_bwa \
    ./"${file}_1_trimmed.fastq.gz" ./"${file}_2_trimmed.fastq.gz" \
    > ../../mapping/"${file}_bwa_dna_pe.sam"

    cd ../../mapping

    # count unique correctly mapped reads
    num=$(samtools view -@ 25 -F 4 -f 2 -h "${file}_bwa_dna_pe.sam"  grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -c)
    echo "${file}_bwa_dna_pe.sam $num" >> stats.txt

    # filter unique correctly mapped reads
    samtools view -@ 25 -F 4 -f 2 -h "${file}_bwa_dna_pe.sam" | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > "${file}_bwa_dna_pe.bam"

    samtools sort -@ 25 "${file}_bwa_dna_pe.bam" -o "${file}_bwa_dna_pe_sorted.bam"
    samtools index -@ 25 "${file}_bwa_dna_pe_sorted.bam"

    echo "finished ${file}"
done