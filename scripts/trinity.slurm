#!/bin/sh
#SBATCH --job-name=trinity
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=10G 
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16


module purge
module load Trinity/2.9.1-foss-2019b-Python-3.7.4

STR1=$1
STR2=$2


#STR="$(ls *.fastq -x1)"

Trinity --seqType fq --left "$1" --right "$2" --max_memory 160G --CPU "$SLURM_CPUS_PER_TASK" --full_cleanup




