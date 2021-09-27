#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=6
##SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
##SBATCH --partition=normal


R1=$1
R2=$2

module purge
module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
trim_galore $R1 $R2 --paired -j 6



