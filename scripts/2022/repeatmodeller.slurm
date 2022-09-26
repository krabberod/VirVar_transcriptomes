#!/bin/sh
#SBATCH --job-name=repeatmodeller
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
##SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
##SBATCH --partition=bigmem

REF=$1
module purge
module load RepeatModeler/2.0.2a-foss-2020b

START=$(date +%s)
echo "###################################################"
echo "# Started: "
echo "#     $(date)" 
echo "###################################################"
echo "# Ref file: $REF"
echo "###################################################"
echo ""

RepeatModeler -database $REF -pa 10 -LTRStruct -rmblast_dir /cluster/software/RMBlast/2.9.0-gompi-2019b/bin 

# Computing runtime
secs=$(($(date +%s)-$START))
echo ""
echo "######################"
echo "Script finished: "
echo "    $(date)"
echo "Running time:"
#echo "#     $(($(date +%si)-$START)) seconds"
printf '%dd:%dh:%02dm:%02ds\n' $(($secs/86400)) $(($secs%86400/3600)) $(($secs%3600/60)) $(($secs%60)) 
echo "######################"

