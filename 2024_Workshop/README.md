# Repository for the 2024 Workshop on VirVar transcriptomes
In our last session, we looked at the virus gene expression. Now it's time to look at the host's response to different kinds of viruses. 

The *Prymnesium kappa* genome is not yet ready, but we have *Haptolina ericina*.

# H.ericina annotation
The H.ericina genome has been annotated with the following tools:
- [Repeatmodeler](https://www.repeatmasker.org/RepeatModeler/): A tool that identifies repetitive elements in eukaryotic genomes. It uses [RECON](https://www.pnas.org/doi/10.1073/pnas.1921046117) and [RepeatScout](https://www.pnas.org/doi/10.1073/pnas.1921046117) to build a repeat library and then uses [Repeatmasker](https://www.pnas.org/doi/10.1073/pnas.1921046117) to mask repetitive elements in the genome. [Flynn et al (2020)](https://www.pnas.org/doi/10.1073/pnas.1921046117)
- Braker2: A pipeline that combines GeneMark-ES and Augustus to predict genes in eukaryotic genomes.
- [Genemark-ES](http://exon.gatech.edu/GeneMark/eukaryotes/): A self-training gene prediction tool that uses unsupervised training to predict genes in eukaryotic genomes.
- [Augustus](http://bioinf.uni-greifswald.de/augustus/): A gene prediction tool that uses a hidden Markov model to predict genes in eukaryotic genomes.

The output from the annotation is available in the folder *Haptolina_ericina_var_UIO028_V1_genome_annotations*, with the following subfolders:
- **01_repeatmasker**: Contains the output from Repeatmasker.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.cat.gz*: Detailed output of the repeatmasker analysis.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.out*: a list of the repetitive elements identified in the genome.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.out.gff*:  a GFF file of the repetitive elements identified in the genome.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.out.html*: a summary of the repetitive elements identified in the genome viewable in a web browser.
  - *Haptolina_ericina_var_UIO028.mainGenome.fasta.tbl*:  Summary statistics of the repetitive elements identified in the genome.
- **02_augustus_proteins**: Contains the predicted proteins from Genemark/Augustus.
  - The gene prediciton was done on the soft masked masked genome.
  - 113,522 predicted genes. Some from mulitple isoforms. 
- **03_augustus_proteins_interproscan**: Contains Interproscan results for the predicted proteins from Augustus. Refer to the [Interproscan documentation](https://www.ebi.ac.uk/interpro/) for more information. Paysan-Lafosse, T., et al. (2022). https://doi.org/10.1093/nar/gkac993. Interproscan is a tool that scans protein sequences against multiple databases to identify conserved domains and functional motifs.
- **04_Eggnog**: Contains the Eggnog results for the predicted proteins from Augustus. Refer to the [Eggnog documentation](http://eggnog5.embl.de/) for more information. Huerta-Cepas, J., et al. (2019). https://doi.org/10.1093/nar/gky1085. Eggnog is a tool that assigns orthology and functional annotations to proteins based on evolutionary relationships.


# Workshop Goals
H.ericina has been infected with two different viruses, and we have the transcriptomes of the host and the viruses for several time points.
1.	PkV-RF01 - H.ericina 028
2.	HeV-RF02 - H.ericina 028

The primary objective is to determine whether the host's response varies when infected with different viruses. Specifically, we want to know if there's a difference in the host's response when H. ericina is infected with either HeV-RF02 or PkV-RF01.

- Genemark-ES was used to predict the genes in the virus HeV-RF02. It was already done for PkV-RF01.
- To map the reads to the host genome, we will use STAR. Previously we used Kallisto. Kallisto and Star aligner are two tools used for RNA-seq read alignment. The main difference between these tools is that Kallisto directly maps the reads to the transcripts, while Star aligner maps the reads to the genome. To use Star aligner, a gtf file is required for genome annotations. One benefit of using Star aligner is that it can map reads to regions of the genome where the annotations have failed. Star is splice aware (in comparison to Bowtie2, which is not).
 
### 1. Create a Combined Reference Genome
First, create a combined reference genome that includes both the host and the viral genomes. This ensures that reads from both organisms can be properly aligned in a single run. This will be done twice, each time with a different virus (but the same host).

Before concatenating the genomes, make sure that the chromosome/contig names are unique between the host and viral genomes. If they overlap, you might need to rename them to avoid conflicts.

Also make sure that the GTF files for the host and viral annotations are formatted correctly. The gene_id and transcript_id attributes should be unique across the combined genome. At this time it is also a good idea to name the genes in the GTF file in a sensible way. I renmaed the genes in the PkV-RF01 GTF file to PkV-RF01_g1, PkV-RF01_g2 etc (see [rename_PkV-RF01_gtf.md](rename_PkV-RF01_gtf.md) for the sed command used).

Concatenate the host and viral genome FASTA files. Using PkVRF01 and H.ericina 028 as an example:
```bash
# cat host_genome.fasta virus_genome.fasta > combined_genome.fasta
cat Haptolina_ericina_var_UIO028.mainGenome.softmasked.fasta PkV-RF01_final.fasta > PkVRF01-He028.fasta
# Similarly for HeVRF02 and H.ericina 028
# cat Haptolina_ericina_var_UIO028.mainGenome.softmasked.fasta HeVRF02.fasta > HeVRF02_He028.genome.fasta
```
Concatenate the GTF files for virus and host annotations:
```bash
# cat host_annotations.gtf virus_annotations.gtf > combined_annotations.gtf
cat Haptolina_ericina_var_UIO028.mainGenome.softmasked.gtf PkV-RF01_final.gtf > PkVRF01-He028.gtf
# Similarly for HeVRF02 and H.ericina 028
```

### 2. Index the Reference Genomes
Index the combined reference genome with STAR. This step only needs to be done once for each combined genome.

STAR is a module on SAGA
```bash
# STAR is a module on SAGA
module load STAR/2.7.11a-GCC-12.3.0
```
Run the following command to index the combined genome:

```bash
STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --runMode genomeGenerate \
     --genomeDir PkVRF01-He028_star \
     --genomeFastaFiles PkVRF01-He028.fasta \
     --sjdbGTFfile PkVRF01-He028.gtff
```
**NB** Set the --sjdbOverhang parameter to 1 less than the read length you are using for the mapping IF the reads are shorter than 100. I tested with 100 bp and the average read length of 136bp. It didn't seem to make a difference.
> The `--sjdbOverhang` is used only at the genome generation step, and tells STAR how many bases to concatenate from donor and acceptor sides of the junctions. If you have 100b reads, the ideal value of `--sjdbOverhang` is 99, which allows the 100bp read to map 99bp on one side, 1bp on the other side. One can think of `--sjdbOverhang` as the maximum possible overhang for your reads.

### 3. Map the Reads

Star requires the genome index to be in a directory with the name of the genome index.  

For quick unzipping of the reads, use the pigz module on SAGA:
```bash
module load pigz/2.7-GCCcore-12.2.0
pigz -d *fq.gz
```


Map the reads to the combined reference genome with STAR:
(if reads are gzipped use **--readFilesCommand zcat**)

```bash
STAR --genomeDir $GENOME_INDEX \
     --readFilesIn read1.fastq read2.fastq \
     --runThreadN $SLURM_CPUS_PER_TASK \
     --outFileNamePrefix outputPrefix \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts
```
Some premlimary results: 
- The % of reads mapping uniquely is a bit low: about 30%. 
- The % of reads mapping to multiple locations is high, higher than "normal" (about 40%).
- Need to check the mapping quality, or other reasons that could explain this.
  - Genome duplication?
  - Polyploidy?

**UPDATE 13 MAY**
The number of reads mapping to multiple locations went down after removing rDNA from the RNA-reads as expected (reads mapped to multiple loci: 3.16%, of reads mapped to too many loci: 1.78%): The percentage of reads mapping uniquely went up a bit to 45.9%, but that is still a bit low. Next step is to adjust parameters that control the minimum number of matched bases to consider a read as mapped and the minimum fraction of the read length that needs to be aligned. I.e. the --outFilterMatchNmin and --outFilterScoreMinOverLread parameters.

> The default behavior for STAR regarding reads that are too short to be mapped is determined by a couple of key parameters:
> 
> --outFilterMatchNmin: This parameter sets the minimum number of matched bases to consider a read as mapped. By default, STAR does not explicitly set this parameter, which means it is essentially governed by the overall alignment criteria set by other parameters.
> 
> --outFilterScoreMinOverLread: This parameter, along with --outFilterMatchNminOverLread, sets the minimum fraction of the read length that needs to be aligned. The default value for --outFilterScoreMinOverLread is 0.66 (or 66%).
> 
> --outFilterMatchNminOverLread: The default for this parameter is 0.66 (or 66%), meaning that at least 66% of the read length needs to match the reference for it to be considered mapped.
> 
> Given these defaults, reads that are shorter than two-thirds of their original length in terms of matched bases will be considered too short to map. This effectively means that if a read is trimmed or inherently short, it might fall below this threshold and be classified as "too short."

Testing with a lower setting for --outFilterMatchNmin and --outFilterScoreMinOverLread. 

```bash
STAR --genomeDir $GENOME \
     --readFilesIn $1 $2 \
     --runThreadN $SLURM_CPUS_PER_TASK \
     --outFileNamePrefix $STR.PkVRF01-He028 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outFilterMatchNminOverLread 0.3 \
     --outFilterScoreMinOverLread 0.3
```
The percentage of reads mapping uniquely went up to 56.88% Which is better (but lowering the mach criteria might also increase the number of false positives).


### 4. Analyze Differential Expression
The resulting **ReadsPerGene.out.tab** will contain counts for both host and viral genes. This data can be loaded into R for differential expression analysis or other downstream analysis. 


### Additional Considerations
- Filtering and Quality Control: Depending on the ratio of host to viral RNA, some alignments might be ambiguous. It might be necessary to use more stringent alignment parameters or post-alignment filtering to distinguish between lowly expressed host transcripts and viral transcripts.
- Viral Transcript Variability: Viruses can have high mutation rates. If applicable, consider allowing more mismatches during the alignment or use tools designed to handle high variability.
- Viral Integration? reads that span host-virus junctions? 
- Chimeric Reads Identification: Look for chimeric reads that span both host and viral sequences. This can be indicative of viral integration into the host genome. Tools like STAR can identify chimeric reads during the alignment process.
- Integration Sites Analysis: Potential tools: VirusSeq (https://pubmed.ncbi.nlm.nih.gov/23162058/), VirusFinder (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064465), or ViFi (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6283451/). These tools can detect viral integration by analyzing chimeric reads and discordant read pairs.
- Host-Virus Interaction: Are there specific pathways or genes that are known to be involved in the host's response to viral infection? 
- The number of genes predicted from the Haptolina ericina genome is quite high. 


