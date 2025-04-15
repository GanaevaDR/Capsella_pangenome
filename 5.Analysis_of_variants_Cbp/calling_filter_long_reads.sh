#!/bin/bash

for i in ./*_sorted.bam; do
    file=$(basename "$i" _sorted.bam)

    bcftools mpileup --threads 30 -f ../reference/Co_Cr_reference.fasta ./"${file}_sorted.bam" | \
    bcftools call -mv -Ob --ploidy 2 \
    -o "${file}_calls.bcf"

    touch "${file}_calls_stats.txt"

    bcftools stats./"${file}_calls.bcf" > ./"${file}_calls_stats.txt"

    bcftools filter -e 'INFO/DP < 5 || QUAL < 30' \
    -Ob ./"${file}_calls.bcf" \
    > ./"${file}_filtered_calls.bcf"

    touch "${file}_filtered_calls_stats.txt"

    bcftools stats ./"${file}_filtered_calls.bcf" > ./"${file}_filtered_calls_stats.txt"

    bcftools index "${file}_calls.bcf"
    bcftools index "${file}_filtered_calls.bcf"

done