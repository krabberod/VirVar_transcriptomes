Contains scripts for transcriptome analysis:
- Slurm scripts for running on cluster
  - *trim_galore.slurm* - for trimming of reads. Trim Galore is a wrapper of Cutadapt with additional triming. https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
  -
- R scripts for analysis

Download the reference genome with annotation from:
```
wget https://raw.githubusercontent.com/RomainBlancMathieu/PkV-RF01/master/PkV-RF01_final.gff
wget https://raw.githubusercontent.com/RomainBlancMathieu/PkV-RF01/master/PkV-RF01_final.fnn
wget https://raw.githubusercontent.com/RomainBlancMathieu/PkV-RF01/master/PkV-RF01_final.fasta
wget https://raw.githubusercontent.com/RomainBlancMathieu/PkV-RF01/master/PkV-RF01_final.faa
```

## STEP 1. Cleaning with trimgalore
For-loop for starting the trim_galore script for all samples.
```
for f in Sample_*; do cp trim_galore.slurm $f; cd $f; sbatch trim_galore.slurm *R1_001.fastq* *R2_001.fastq*; cd ..; done
```
Copy the trimmed reads to a new folder and keep the folder structure:
```
rsync -avR 00_data/Sample_*/*val* trimmed
```

This loop creates new directories and softlinks to the trimmed reads:
```
for f in Sample_*; do
mkdir -p ./../06_genome_mapping/$f;
REF=$(readlink -f $f);
echo $REF;
BACK=$(pwd);
cd ./../06_genome_mapping/$f;
ln -s $REF/*val* .;
cd $BACK;
done
```

## STEP 2. Mapping with bowtie2 or kallisto
Several mapping programs can be used, but for now we will use bowtie to map to the genome and kallisto to map to the orf (or predicted genes). Kallisto allows quantification.

### 2a Mapping to the genome with Bowtie2
Mapping requires that the reference is indexed. i.e.
```
module purge
ml Bowtie2/2.4.4-GCC-10.3.0
bowtie2-build PkV-RF01_final.fasta PkV-RF01_final_genome
```
Set the reference in the bowtie.slurm script and run it for each sample with a for-loop

### 2b Mapping to predicted genes with Kallisto
Mapping the reads to the CDS using Kallisto.

https://pachterlab.github.io/kallisto/manual
https://pachterlab.github.io/kallisto/starting
```
module purge
ml kallisto/0.46.1-foss-2020a
kallisto index PkV-RF01_final.fnn -i PkV-RF01_final_cds
```
Set the reference in the kallisto.slurm script and run it for each sample


### Step 4. R analysis
See the script in the R_data folder. It's based on the DESeq2 package in R. 
http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
