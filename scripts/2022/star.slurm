#!/bin/sh
#SBATCH --job-name=star
#SBATCH --account=nn9699k
#SBATCH --output=slurm-%j.base
##SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=6G
##SBATCH --partition=bigmem


START=$(date +%s)
echo "###################################################"
echo "# Started: "
echo "#     $(date)" 
echo "###################################################"
echo ""
module purge
module load STAR/2.7.10a_alpha_220818-GCC-10.3.0 

GENOME=/cluster/work/users/anderkkr/90_Virvar/03_meta_transcriptome/3_PkVRF01_Pk033/00_genome/PkV-RF01_stardb
# /cluster/work/users/anderkkr/90_Virvar/03_meta_transcriptome/1_PkVRF01_He028/04_map_to_genome_pr_sample/HeUIO028_genome_STAR_v2.10a
echo "GENOME is $GENOME"

R1=$1
R2=$2
echo "R1 and R2 are: " $R1 $R2


STAR --genomeDir $GENOME --readFilesIn $R1 $R2 --outSAMtype BAM Unsorted --runThreadN $SLURM_CPUS_PER_TASK --readFilesCommand zcat --outReadsUnmapped Fastx



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

