Contains scripts for transcriptome analysis: 
- Slurm scripts for running on cluster
- R scripts for analysis 

For loop for starting several scripts: 

for f in Sample_*; do cp trim_galore.slurm $f; cd $f; sbatch trim_galore.slurm *R1_001.fastq* *R2_001.fastq*; cd ..; done
