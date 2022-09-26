#!/bin/sh
#SBATCH --job-name=spades
#SBATCH --account=nn9699k
#SBATCH --output=slurm-%j.base
#SBATCH --time=24:0:0
#SBATCH --mem-per-cpu=4G
##SBATCH --partition=bigmem
#SBATCH --cpus-per-task=20


module purge

START=$(date +%s)
echo "###################################################"
echo "# Started: "
echo "#     $(date)"
echo "###################################################"

module load  SPAdes/3.15.3-GCC-10.3.0

STR1=$1
STR2=$2



spades.py -1 $STR1 -2 $STR2 --rna -o rnaspades 



