# The rDNA in He028
- The ribosomal rDNA gene cluster is a highly conserved region in eukaryotic genomes. It is a multi-copy region that is often difficult to assemble and map reads to.
- The rDNA genes have been annotated by an HMM scanning approach. 
- 61 copies were found in the He028 genome, on 13 scaffolds.
- Blasting against PR2 (and Silva) showed possible contamination with a chlorophyte, and/or from chloroplasts.

### Mapping reads to the rDNA
To check for the amount of rDNA reads in the He028 RNAseq data, I mapped the reads to the rDNA genes using bowtie2. 
```bash
module load Bowtie2/2.5.1-GCC-12.2.0
bowtie2-build rDNA.fasta rDNA
DB=/cluster/work/users/anderkkr/140_virvar/01_genomes/He028/ribop/rDNA
bowtie2 -x $DB -1 26-a07va_S190_L004_R1_001_val_1.fq -2 26-a07va_S190_L004_R2_001_val_2.fq -S He028_rDNA.sam --un-conc non_rrna_reads.fq
```