#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=0:10:00
#SBATCH --cpus-per-task=8
##SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=normal
##SBATCH --partition=optimist
##SBACTH --qos=optimist


#STR="$(pwd | sed 's/.*\///g')"

#cp *_R1_* "$STR"_R1.fq.gz;
#cp *_R2_* "$STR"_R2.fq.gz;


module purge
#module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
module load cutadapt/2.10-GCCcore-9.3.0-Python-3.8.2
#trim_galore *_R1_* *_R2* --paired

STR=$(pwd | sed 's/^.*\///')

cutadapt --cores "$SLURM_CPUS_PER_TASK" -a ^CCAGCASCYGCGGTAATTCC...TCATYRATCAAGAACGAAAGT -A ^ACTTTCGTTCTTGATYRATGA...GGAATTACCGCRGSTGCTGG --discard-untrimmed -o "$STR"_cutadapt_R1.fastq.gz -p "$STR"_cutadapt_R2.fastq.gz *R1* *R2*



