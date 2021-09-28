#!/bin/bash
#SBATCH --job-name=bowtie
#SBATCH --account=nn9725k
#SBATCH --output=slurm-%j.base
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=normal

module purge
#module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
ml Bowtie2/2.4.4-GCC-10.3.0
ml SAMtools/1.12-GCC-10.3.0

REF=/cluster/work/users/anderkkr/57_Virvar/00_data/Fasta_genomes/PkV-RF01_genebank

bowtie2 -p 16 --no-unal -x $REF -1 *R1* -2 *R2* > mapped.sam

samtools sort -@ 16 mapped.sam > sorted.bam
samtools index sorted.bam
samtools view -h sorted.bam | cut -f 3,3 | grep 'PkV-RF01' | sort | uniq -c > counts2.txt
