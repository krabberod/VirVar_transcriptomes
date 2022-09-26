
#!/bin/bash
#SBATCH --job-name=Braker2
#SBATCH --output=slurm-%j.base
#SBATCH --account=nn9525k
#SBATCH --time=28:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G


module --quiet purge  
module load Anaconda3/2019.03

source activate /cluster/projects/nn9525k/Programs/conda_envs/braker2

# Export the executables: 
export BSSMDIR=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export GTHDATADIR=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export ALIGNMENT_TOOL_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export BAMTOOLS_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export SAMTOOLS_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export DIAMOND_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export CDBTOOLS_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export PYTHON3_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export AUGUSTUS_BIN_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/bin/
export AUGUSTUS_SCRIPTS_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/git_Augustus/Augustus/scripts/
export MAKEHUB_PATH=/cluster/projects/nn9525k/Programs/conda_envs/braker2/git_Augustus/MakeHub/
export GENEMARK_PATH=/cluster/projects/nn9525k/Programs/gmes_linux_64/
# Remember to get key to GeneMark!
# Config_path is writable:

export AUGUSTUS_CONFIG_PATH=/cluster/work/users/anderkkr/90_Virvar/02_genome_annotation/01_braker/augustus_config
echo $AUGUSTUS_CONFIG_PATH
# Copy all config-files from augustus to your folders, and make it writeable:  If this is not done already
# cp -r /cluster/software/AUGUSTUS/3.4.0-foss-2020b/config/ $AUGUSTUS_CONFIG_PATH
# rnaseq.bam is the mapped reads to the genome. 

braker.pl --genome=Haptolina_ericina_var_UIO028.mainGenome.fasta --bam=rnaseq.bam --cores $SLURM_CPUS_PER_TASK --workingdir=braker_prot --makehub --email a.k.krabberod@ibv.uio.no
