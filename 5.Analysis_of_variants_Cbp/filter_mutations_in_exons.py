import pandas as pd

def read_gtf(gtf_file):
    """Read the GTF file and extract exon coordinates."""
    exon_coords = []
    
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue  
            fields = line.strip().split('\t')
            if fields[2] == 'exon':  # Check if the feature is an exon
                chromosome = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                exon_coords.append((chromosome, start, end))
    
    return exon_coords

def filter_vcf(vcf_file, exon_coords, vcf_out):
    """Filter the VCF file based on exon coordinates."""
    exon_set = set()
    
    for chrom, start, end in exon_coords:
        exon_set.add((chrom, start, end))
    
    with open(vcf_file, 'r') as vcf, open(vcf_out, 'w') as output_vcf:
        for line in vcf:
            if line.startswith('#'):
                output_vcf.write(line)  
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            
            # Check if the position is within any of the exon ranges
            for (ex_chrom, ex_start, ex_end) in exon_set:
                if chrom == ex_chrom and ex_start <= pos <= ex_end:
                    output_vcf.write(line)
                    break  


# Provide files and run the functions
gtf_file = 'Co_Cr_reference.gtf' 
vcf_file = 'Eng_60_calls_filt_5_anno_target_high_damaging.ann.vcf'   
vcf_out = 'Eng_exons_damaged.vcf'
exon_coordinates = read_gtf(gtf_file)
filter_vcf(vcf_file, exon_coordinates, vcf_out)

gtf_file = 'Co_Cr_reference.gtf' 
vcf_file = 'KBG_60_calls_filt_5_anno_target_high_damaging.ann.vcf'   
vcf_out = 'KBG_exons_damaged.vcf'
filter_vcf(vcf_file, exon_coordinates, vcf_out)

gtf_file = 'Co_Cr_reference.gtf'  
vcf_file = 'MSK_60_calls_filt_5_anno_target_high_damaging.ann.vcf'   
vcf_out = 'MSK_exons_damaged.vcf'
filter_vcf(vcf_file, exon_coordinates, vcf_out)

gtf_file = 'Co_Cr_reference.gtf'  
vcf_file = 'Mur_60_calls_filt_5_anno_target_high_damaging.ann.vcf'   
vcf_out = 'Mur_exons_damaged.vcf'
filter_vcf(vcf_file, exon_coordinates, vcf_out)

gtf_file = 'Co_Cr_reference.gtf' 
vcf_file = 'le3_60_calls_filt_5_anno_target_high_damaging.ann.vcf'   
vcf_out = 'le3_exons_damaged.vcf'
filter_vcf(vcf_file, exon_coordinates, vcf_out)