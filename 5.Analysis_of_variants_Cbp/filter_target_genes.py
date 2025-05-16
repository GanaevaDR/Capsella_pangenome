#!/usr/bin/env python3

import re

def filter_vcf_by_genes(genes_file, vcf_file, output_file):
    """
    Filters a VCF file to include only lines that contain specified gene IDs.

    Parameters:
    - genes_file: str, path to the file containing gene IDs of interest.
    - vcf_file: str, path to the input VCF file to be filtered.
    - output_file: str, path to the output file where filtered results will be saved.
    """

    # Read gene IDs into a set for fast lookup
    with open(genes_file, 'r') as f:
        genes_of_interest = set(line.strip() for line in f)

    # Compile a regex pattern for filtering
    pattern = re.compile(r'GENE=({})'.format('|'.join(re.escape(gene) for gene in genes_of_interest)))

    # Filter the VCF file
    with open(vcf_file, 'r') as vcf, open(output_file, 'w') as out:
        for line in vcf:
            if pattern.search(line):
                out.write(line)

    print(f"Filtering complete for {vcf_file}. Results saved to {output_file}.")

# Example use
genes_file = 'Conserved_genes.txt'
vcf_file = "Eng_60_calls_anno.vcf"
output_file = 'Eng_60_calls_filt_5_anno_target_genes.vcf'

filter_vcf_by_genes(genes_file, vcf_file, output_file)