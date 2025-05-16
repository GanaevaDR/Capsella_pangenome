# Capsella_pangenome

This project is dedicated to the analysis of gene evolution of recent allopolyploid plant _Capsella bursa-pastoris_.


## 1. Orthologization
- ProteinOrtho (by default)
  
``` proteinortho --cpus=1 C_rub_proteins.faa C_ori_proteins.faa > proteinortho.log ```

- ProteinOrtho, synteny mode
  
``` proteinortho --synteny --cpus=1 C_rub_proteins.faa C_ori_proteins.faa > proteinortho_synteny.log ```

- OrthoFinder
  
``` orthofinder.py -f ./ -t 10 -o ./orthofinder_res ```

- SynGAP
```
syngap dual \
--sp1fa=C_orientalis.fasta \
--sp1gff=C_ori_longest_transcripts_cor_id.gff \
--sp2fa=C_rubella.fasta \
--sp2gff=C_rub_longest_transcripts_cor_id.gff \
--sp1=Co \
--sp2=Cr \
--datatype prot
```

```
 syngap genepair \
--sp1fa=Co.fa \
--sp1gff=Co.SynGAP.gff3 \
--sp2fa=Cr.fa \
--sp2gff=Cr.SynGAP.gff3 \
--sp1=Co \
--sp2=Cr \
--datatype prot
```

- **extract_orthopairs_syngap.R** - parse output of SynGAP to extract oprthopairs
- **cds_length_estimation.R** - calculate CDS length difference in orthopairs
- **cds_length_plots.R** - draw CDS length distribution plots

- minimap2 commands for alignment of CDS to the genome:

```
#create index files:
minimap2 -d C_orientalis.mmi C_orientalis.fasta
minimap2 -d C_rubella.mmi C_rubella.fasta
```
```
# run minimap2 with permission for intron gaps:
minimap2 -ax splice -k14 -uf -t 10 -G 10k C_orientalis.fasta Cr_larger.fasta > Cr_larger_genes_to_Co_genome.sam
minimap2 -ax splice -k14 -uf -t 10 -G 10k C_rubella.fasta Co_larger.fasta > Co_larger_genes_to_Cr_genome.sam
```

## 2. Processing of DNA sequencing data
- **raw_reads_qc.sh** - perform quality controls of raw DNA-seq reads
- **trim_qc_short_reads.sh** - perform trimming and quality control of trimmed reads
- **trimming_long_reads.sh** - trimming of long reads
- **mapping_filter_short_reads.sh** - mapping and filtration of unique correctly mapped Illumina short reads
- **mapping_filter_long_reads.sh** - mapping and filtration of PacBio long reads
- **gene_coverage.sh** - calculate average gene coverage (depth) per sample
- **gene_breadth.sh** - calculate average gene breadth per sample
  
## 3. Analysis of conserved genes in parent species 
**conserved_genes_parent_species.R** 

This script:

- merges files with average gene coverage per sample into a single table
- filters out low covered samples
- filters out genes with low coverage
- outputs conserved genes that belong to orthopairs in paren species

## 4. Analysis of absent and core genes in _Capsella bursa-pastoris_
**absent_core_genes_Cbp.R**

This script:

- merges files with average gene breadth per sample into a single table
- plots average gene breadth distribution subgenome-wise
- filters genes with average breadth 0 and 1
- provides intersection of genes with average breadth 0 and 1 in samples and draws UpSet diagrams
- filters 1/0 fraction (one gene from the orthopair is present)

## 5. Analysis of nucleotide variants in core genes of _Capsella bursa-pastoris_
- **calling_filter_long_reads.sh** - perform calling and variant filtration for PacBio reads of _Capsella bursa-pastoris_
- **annotate_variants.sh** - merge variants and gene annotation
- **filter_target_genes.py** - filter annotated variants for those belonging to target genes of our interest
- annotate variants using SnpEff:
```
# Build database in SNPEff
snpEff build -gff3 -v Co_Cr_reference
```
**snpeff_run.sh** - run SnpEff on files with variants in VCF format
- **filter_mutations_in_exons.py** - filter records in VCF files that corresponds to mutations in exons
- **filter_unique_conserved_genes.py** - extract unique genes that correspond to conserved parent gene set from VCF files
- **GO_analysis_mutations.R** - perfrom analysis of GO enrichment using gene sets obtained from filter_unique_conserved_genes.py and output from David tool (https://davidbioinformatics.nih.gov/summary.jsp)

## 6. RNA-seq analysis
- **RNA_seq_raw_qc.sh** - quality controls of raw RNA-seq reads
- **RNA_seq_trimming.sh** - perform trimming and quality control of trimmed reads
- **RNA_seq_mapping_quantification.sh** - mapping and quantification of trimmed reads
- **RNA_seq_merge_readcount.py** - merge output tables with readcounts per sample into a single table
