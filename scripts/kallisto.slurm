#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=normal

module purge
ml kallisto/0.46.1-foss-2020a

REF=/cluster/work/users/anderkkr/57_Virvar/00_data/Genome_annotated/PkV-RF01_final_cds
OUT=kallisto_results
R1=*R1_001_val_1*
R2=*R2_001_val_2*


kallisto quant -i $REF -t "$SLURM_CPUS_PER_TASK" -o $OUT $R1 $R2
