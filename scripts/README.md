Contains scripts for transcriptome analysis:
- Slurm scripts for running on cluster
  - *trim_galore.slurm* - for trimming of reads. Trim Galore is a wrapper of Cutadapt with additional triming. https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
  -
- R scripts for analysis

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
mkdir -p ./../04_orf_mapping/$f;
REF=$(readlink -f $f);
echo $REF;
BACK=$(pwd);
cd ./../04_orf_mapping/$f;
ln -s $REF/*val* .;
cd $BACK;
done
```

## STEP 2. Mapping with bowtie2 or kallisto
Several mapping programs can be used, but for now we will use bowtie to map to the genome and kallisto to map to the orf (or predicted genes). Kallisto allows quantification.

### 2a Mapping to the genome with Bowtie2
Mapping requires that the reference is indexed. i.e.
```
ml
bowtie2-build PkV-RF01_genebank.fasta PkV-RF01_genebank
```
Set the reference in the bowtie.slurm script and run it for each sample with a for-loop

### 2b Mapping to predicted genes with Kallisto
```
ml kallisto/0.46.1-foss-2020a
kallisto index PkV-RF01_genebank.fasta -i PkV-RF01_genebank
```
Set the reference in the kallisto.slurm script and run it for each sample
