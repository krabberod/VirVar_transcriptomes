#!/bin/sh
#SBATCH --job-name=cd-hit
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
##SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=108:00:00
#SBATCH --mem-per-cpu=6G
##SBATCH --partition=bigmem


REF=$1
STR=$(echo $REF| sed 's/.fasta//')

START=$(date +%s)
echo "###################################################"
echo "# Started: "
echo "#     $(date)" 
echo "# Inputfile: $REF "
echo "###################################################"
echo ""

module purge
ml CD-HIT/4.8.1-GCC-10.2.0
cd-hit-est -M 0 -i "$REF" -o "$STR"_0.97 -c 0.97 -T 0

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

