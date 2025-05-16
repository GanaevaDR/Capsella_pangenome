#!/bin/bash

for i in ./*_target_genes.vcf; do
    file=$(basename "$i" .vcf)

    snpEff ann -s "${file}_stats.html" Co_Cr_reference "${file}.vcf" > "${file}.ann.vcf"

done