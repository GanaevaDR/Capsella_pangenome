#!/bin/bash
for fold in ./trimmed_data/*_trimmed.fastq.gz; do

    file=$(basename "$fold" _trimmed.fastq.gz)

    STAR --sjdbGTFfile ./reference/Co_Cr_reference.gtf \
    --sjdbGTFfeatureExon exon \
    --sjdbGTFtagExonParentTranscript transcript_id \
    --sjdbGTFtagExonParentGene gene_id \
    --runThreadN 20 --quantMode GeneCounts \
    --genomeDir ./reference/star_index \
    --readFilesIn ./trimmed_data/"${file}_trimmed.fastq.gz" \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix ./mapping/"${file}"
        
done
