#!/bin/bash
#SBATCH --job-name=augustus
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=8
##SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=bigmem

STR=$1
#"$(ls *.fasta -x1)"
module purge
module load AUGUSTUS/3.3.2-intel-2018b-Python-2.7.15

#augustus --species=arabidopsis "$STR" > "$STR".augustus
augustus --species=chicken "$STR" > "$STR".augustus
getAnnoFasta.pl "$STR".augustus --seqfile="$STR"



