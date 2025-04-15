#!/bin/bash

for i in ./*.fastq.gz; do

    file=$(basename "$i" .fastq.gz)
    cutadapt --poly-a -o ./"${file}_trimmed.fastq.gz"  ./"${file}.fastq.gz"

    echo "finished ${file}"
    
done