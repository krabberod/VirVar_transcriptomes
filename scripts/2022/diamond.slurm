#!/bin/sh
#SBATCH --job-name=diamond
#SBATCH --account=nn9244k
#SBATCH --output=slurm-%j.base
#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=6G

module -q purge
ml DIAMOND/2.0.13-GCC-11.2.0

FASTA=$1
# DB=/cluster/projects/nn9725k/Programmer/databases/tara_SMAGs/proteins/SMAGs_prot.dmnd
# DB=/cluster/projects/nn9725k/Programmer/databases/uniprot_swissprot/uniprot.dmnd
# DB=/cluster/projects/nn9725k/Programmer/databases/refseq_plants/refseq.plants.dmnd
DB=/cluster/work/users/anderkkr/tmp/nr.dmnd

STR=$(basename -s .fasta $FASTA)
STR2=$(basename $DB)
OUT_TSV="$STR"."$STR2".diamond.tsv

START=$(date +%s)
echo "###################################################"
echo "# Started: "
echo "#     $(date)"
echo "###################################################"
echo "# REF: $1 "
echo ""


diamond blastp -p $SLURM_CPUS_PER_TASK --db $DB -q $FASTA --out $OUT_TSV

# Computing runtime
secs=$(($(date +%s)-$START))
echo ""
echo "######################"
echo "Script finished: "
echo "    $(date)"
echo "Running time:"
printf '%dd:%dh:%02dm:%02ds\n' $(($secs/86400)) $(($secs%86400/3600)) $(($secs%3600/60)) $(($secs%60))
echo "######################"
