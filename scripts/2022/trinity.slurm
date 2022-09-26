#!/bin/sh
#SBATCH --job-name=trinity
#SBATCH --account=nn9699k
#SBATCH --output=slurm-%j.base
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=5G 
##SBATCH --partition=bigmem
#SBATCH --cpus-per-task=30


module purge
# module load Trinity/2.9.1-foss-2019b-Python-3.7.4
module load Trinity/2.10.0-foss-2020a-Python-3.8.2

START=$(date +%s)
echo "###################################################"
echo "# Started: "
echo "#     $(date)"
echo "###################################################"
echo "# R1: $1"
echo "# R2: $2"
echo "###################################################"
echo ""


#STR="$(ls *.fastq -x1)"

Trinity --seqType fq --left "$1" --right "$2" --max_memory 150G --CPU "$SLURM_CPUS_PER_TASK" --full_cleanup 


# Computing runtime
secs=$(($(date +%s)-$START))
echo ""
echo "######################"
echo "Script finished: "
echo "    $(date)"
echo "Running time:"
printf '%dd:%dh:%02dm:%02ds\n' $(($secs/86400)) $(($secs%86400/3600)) $(($secs%3600/60)) $(($secs%60))
echo "######################"

