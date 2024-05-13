#!/bin/sh
#SBATCH --job-name=star_index
#SBATCH --account=nn9845k
#SBATCH --output=slurm-%j.base
##SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=6G
##SBATCH --partition=bigmem

set -o errexit
set -o nounset

START=$(date +%s)
echo "###################################################"
echo "# Started: "
echo "#     $(date)" 
echo "###################################################"
echo ""

#SLURM_CPUS_PER_TASK=2
module purge
module load STAR/2.7.11a-GCC-12.3.0
STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --runMode genomeGenerate \
     --genomeDir PkVRF01-He028_star \
     --genomeFastaFiles PkVRF01-He028.fasta \
     --sjdbGTFfile PkVRF01-He028.gtf \
     --genomeSAindexNbases 12

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

