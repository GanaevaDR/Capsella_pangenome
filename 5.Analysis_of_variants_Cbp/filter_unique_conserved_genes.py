#!/usr/bin/env python3

import re

def filter_vcf_by_genes(genes_file, vcf_file, output_genes_file):
    """
    Filters a VCF file to include only lines that contain specified gene IDs,
    and outputs a list of unique gene IDs found in the VCF file.

    Parameters:
    - genes_file: str, path to the file containing gene IDs of interest
    - vcf_file: str, path to the input VCF file to be filtered
    - output_genes_file: str, path to the output file where found gene IDs will be saved
    """

    # Read gene IDs into a set
    with open(genes_file, 'r') as f:
        genes_of_interest = set(line.strip() for line in f)

    # Compile a regex pattern for filtering
    pattern = re.compile(r'GENE=({})'.format('|'.join(re.escape(gene) for gene in genes_of_interest)))

    # Set to store found genes
    found_genes = set()

    # Filter the VCF file
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if pattern.search(line):
                # Extract gene names from the line and add to found_genes
                matches = pattern.findall(line)
                found_genes.update(matches)

    # Write found genes to the output genes file
    with open(output_genes_file, 'w') as genes_out:
        for gene in sorted(found_genes):
            genes_out.write(f"{gene}\n")

    print(f"Filtering complete for {vcf_file}.")
    print(f"Found genes saved to {output_genes_file}.")

# File paths
genes_file = 'conserved_genes.txt'
vcf_file = "Eng_exons_damaged.vcf"
output_genes_file = "Eng_damaged_exons.txt"
filter_vcf_by_genes(genes_file, vcf_file, output_genes_file)


genes_file = 'conserved_genes.txt'
vcf_file = "KBG_exons_damaged.vcf"
output_genes_file = "KBG_damaged_exons.txt"
filter_vcf_by_genes(genes_file, vcf_file, output_genes_file)


genes_file = 'conserved_genes.txt'
vcf_file = "MSK_exons_damaged.vcf"
output_genes_file = "MSK_damaged_exons.txt"
filter_vcf_by_genes(genes_file, vcf_file, output_genes_file)


genes_file = 'conserved_genes.txt'
vcf_file = "Mur_exons_damaged.vcf"
output_genes_file = "Mur_damaged_exons.txt"
filter_vcf_by_genes(genes_file, vcf_file, output_genes_file)


genes_file = 'conserved_genes.txt'
vcf_file = "le3_exons_damaged.vcf"
output_genes_file = "le3_damaged_exons.txt"
filter_vcf_by_genes(genes_file, vcf_file, output_genes_file)


