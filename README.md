# Capsella_pangenome

This project is dedicated to the analysis of gene evolution of recent allopolyploid plant _Capsella bursa-pastoris_


### Orthologization
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

``` syngap genepair \
--sp1fa=Co.fa \
--sp1gff=Co.SynGAP.gff3 \
--sp2fa=Cr.fa \
--sp2gff=Cr.SynGAP.gff3 \
--sp1=Co \
--sp2=Cr \
--datatype prot ```


### Processing of DNA sequencing data


### Analysis of conserved genes in parent species 


### Analysis of absent and core genes in _Capsella bursa-pastoris_


### Analysis of nucleotide variants in core genes of _Capsella bursa-pastoris_
