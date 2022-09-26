#!/bin/sh
#SBATCH --job-name=fastp
#SBATCH --account=nn9699k
#SBATCH --output=slurm-%j.base
##SBATCH --nodes=1 
#SBATCH --cpus-per-task=8
#SBATCH --time=1:30:00
#SBATCH --mem-per-cpu=10G
##SBATCH --partition=bigmem


module purge
ml fastp/0.23.2-GCC-11.2.0

STR1=$1
STR2=$2

# Output
STR3=$(echo $1 | sed 's/fq.gz/trim_fq.gz/')
STR4=$(echo $2 | sed 's/fq.gz/trim.fq.gz/')

fastp -i $STR1 -I $STR2 -o $STR3 -O $STR4 -w "$SLURM_CPUS_PER_TASK" -l 100
