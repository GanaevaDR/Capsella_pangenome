#!/bin/bash

for i in ./*_1.fastq.gz; do

        file=$(basename "$i" _1.fastq.gz)
        fastp --in1  ./"${file}_1.fastq.gz"  --in2 ./"${file}_2.fastq.gz" \
        --out1 ./"${file}_1_trimmed.fastq.gz" --out2 ./"${file}_2_trimmed.fastq.gz"  \
        -l 50 --cut_front --cut_right --trim_poly_x --dedup \
        --adapter_fasta adapters.all  \
        --thread 6 -h ./"${file}_trim.html"

        echo "finished ${file}"

done

fastqc -o ./fastqc_trimmed/  ./*_trimmed.fastq.gz
multiqc  . -o  ./fastqc_trimmed