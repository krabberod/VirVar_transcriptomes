# Repository for the 2024 Workshop on VirVar transcriptomes
In our last session, we looked at the virus gene expression. Now it's time to look at the host. 

The *Prymnesium kappa* genome is not yet ready, but we have *Haptolina ericina*.

### H.ericina annotation
The H.ericina genome has been annotated with the following tools:
- [Repeatmodeler](https://www.repeatmasker.org/RepeatModeler/): A tool that identifies repetitive elements in eukaryotic genomes. It uses [RECON](https://www.pnas.org/doi/10.1073/pnas.1921046117) and [RepeatScout](https://www.pnas.org/doi/10.1073/pnas.1921046117) to build a repeat library and then uses [Repeatmasker](https://www.pnas.org/doi/10.1073/pnas.1921046117) to mask repetitive elements in the genome. [Flynn et al (2020)](https://www.pnas.org/doi/10.1073/pnas.1921046117)
- Braker2: A pipeline that combines GeneMark-ES and Augustus to predict genes in eukaryotic genomes.
- [Genemark-ES](http://exon.gatech.edu/GeneMark/eukaryotes/): A self-training gene prediction tool that uses unsupervised training to predict genes in eukaryotic genomes.
- [Augustus](http://bioinf.uni-greifswald.de/augustus/): A gene prediction tool that uses a hidden Markov model to predict genes in eukaryotic genomes.

The output from the annotation is available in the folder *Haptolina_ericina_var_UIO028_V1_genome_annotations*, with the following subfolders:
- **01_repeatmasker**: Contains the output from Repeatmasker.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.cat.gz* : Detailed output of the repeatmasker analysis.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.out* : a list of the repetitive elements identified in the genome.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.out.gff* :  a GFF file of the repetitive elements identified in the genome.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.out.html* : a summary of the repetitive elements identified in the genome viewable in a web browser.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.tbl*:  Summary statistics of the repetitive elements identified in the genome.
- **02_augustus_proteins**: Contains the predicted proteins from Genemark/Augustus.
  - 113,522 predicted genes. Some from mulitple isoforms. 
- **03_augustus_proteins_interproscan**: Contains Interproscan results for the predicted proteins from Augustus. Refer to the [Interproscan documentation](https://www.ebi.ac.uk/interpro/) for more information. Paysan-Lafosse, T., et al. (2022). https://doi.org/10.1093/nar/gkac993. Interproscan is a tool that scans protein sequences against multiple databases to identify conserved domains and functional motifs.
- **04_Eggnog**: Contains the Eggnog results for the predicted proteins from Augustus. Refer to the [Eggnog documentation](http://eggnog5.embl.de/) for more information. Huerta-Cepas, J., et al. (2019). https://doi.org/10.1093/nar/gky1085. Eggnog is a tool that assigns orthology and functional annotations to proteins based on evolutionary relationships.


# Workshop Goals
H.ericina has been infected with two different viruses, and we have the transcriptomes of the host and the viruses for several time points.
1.	Pk-VRF01 - H.ericina 028
2.	He-VRF02 - H.ericina 028

The primary objective is to determine whether the host's response varies when infected with different viruses. Specifically, we want to know if there's a difference in the host's response when H. ericina is infected with either HeVRF02 or PkVRF01.

- Genemark-ES was used to predict the genes in the virus He-VRF02. It was already done for Pk-VRF01.
- To map the reads to the host genome, we will use STAR.
  
### 1. Create a Combined Reference Genome
First, create a combined reference genome that includes both the host and the viral genomes. This ensures that reads from both organisms can be properly aligned in a single run. This will be done twice, each time with a different virus.

Before concatenating the genomes, make sure that the chromosome/contig names are unique between the host and viral genomes. If they overlap, you might need to rename them to avoid conflicts.

Concatenate the host and viral genome FASTA files:
```bash
# cat host_genome.fasta virus_genome.fasta > combined_genome.fasta
cat PkVR01.fasta H.ericina_028.fasta > PkVRF01_He028.genome.fasta
cat HeVRF02.fasta H.ericina_028.fasta > HeVRF02_He028.genome.fasta
```
Concatenate the GTF files for virus and host annotations:
```bash
# cat host_annotations.gtf virus_annotations.gtf > combined_annotations.gtf
cat PkVR01.gtf H.ericina_028.gtf > PkVRF01_He028.gtf
cat PkVR01.gtf H.ericina_028.gtf > HeVRF02_He028.gtf
```

### 2. Index the Reference Genomes
Index the combined reference genome with STAR:
```bash
STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --runMode genomeGenerate \
     --genomeDir /path/to/combinedGenomeIndex \
     --genomeFastaFiles /path/to/combined_genome.fasta \
     --sjdbGTFfile /path/to/combined_annotations.gtf \
     --sjdbOverhang ReadLength-1
```
### 3. Map the Reads
Map the reads to the combined reference genome with STAR:
```bash
STAR --genomeDir /path/to/combinedGenomeIndex \
     --readFilesIn /path/to/read1.fastq /path/to/read2.fastq \
     --runThreadN $SLURM_CPUS_PER_TASK \
     --outFileNamePrefix /path/to/outputPrefix \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts
```
### 4. Analyze Differential Expression
The resulting **ReadsPerGene.out.tab** will contain counts for both host and viral genes. This data can be loaded into R for differential expression analysis or other downstream analysis. 
```R
library(DESeq2)
countData <- read.delim("/path/to/ReadsPerGene.out.tab", row.names=1)
colData <- data.frame(condition=factor(c("control", "treatment")))
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design= ~ condition)
dds <- DESeq(dds)
res <- results(dds)
```


### Additional Considerations
- Quantifying Viral Integration or Replication: If you're interested in quantifying aspects like viral integration or replication levels, you might need to look specifically at reads that span host-virus junctions or are uniquely mapping to viral sequences.
- Filtering and Quality Control: Depending on the ratio of host to viral RNA, some alignments might be ambiguous. It might be necessary to use more stringent alignment parameters or post-alignment filtering to distinguish between lowly expressed host transcripts and viral transcripts.
- Viral Transcript Variability: Viruses can have high mutation rates. If applicable, consider allowing more mismatches during the alignment or use tools designed to handle high variability.