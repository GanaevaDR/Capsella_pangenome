#!/bin/bash

fastqc -o ./raw_fastqc/  ./raw_data/*.fastq.gz
multiqc  ./raw_data -o  ./raw_fastqc
