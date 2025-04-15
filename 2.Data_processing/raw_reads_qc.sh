#!/bin/bash

mkdir -p fastqc

for project in */; do
        folder=$(basename "$project")
        mkdir ./fastqc/"$folder"
        fastqc -o "./fastqc/$folder/"  "./raw_data/$folder/"*
        multiqc  "./fastqc/$folder" -o  "./fastqc/$folder"
        echo "finished ${folder}"
done