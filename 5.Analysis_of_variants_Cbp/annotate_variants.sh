#!/bin/bash

for i in ./*filtered_calls.bcf; do
    file=$(basename "$i" .bcf)

bcftools view -O v "${file}.bcf" -o "${file}.vcf"

bcftools annotate \
  -a Co_Cr_chr_genes_filt.bed.gz \
  -c CHROM,FROM,TO,GENE \
  -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
  "${file}.vcf"  > "${file}_anno.vcf"

done