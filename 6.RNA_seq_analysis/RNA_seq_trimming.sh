#!/bin/bash

for i in ./*.fastq.gz; do

	file=$(basename "$i" .fastq.gz)
        fastp --in1  ./"${file}.fastq.gz" \
        --out1 ./trimmed_data/"${file}_trimmed.fastq.gz"  \
	--cut_front --cut_right \
        -l 35 --trim_poly_x \
        --adapter_fasta ./adapters/adapters.all  \
        --thread 20 -h ./trimmed_data/"${file}_trim.html"

        echo "finished ${file}"
done
cd ..

fastqc -o ./trimmed_fastqc/  ./trimmed_data/*_trimmed.fastq.gz
multiqc  ./trimmed_fastqc -o  ./trimmed_fastqc