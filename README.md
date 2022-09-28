# VirVar transcriptomes
Resources for Virar. Workshop 27-29 Sept 2021


From Jan:
- FASTQC for visualizing read quality, adaptor presence etc.
- Cutadapt for adapter trimming (you can also use Trimmomatic)
- Trinity for reference transcriptome assembly (all treatments together if you lack a better reference like a genome)
- Bowtie2 for mapping reads to my contig(s) of interest: I used this for a focused analysis of a small subset of molecules
- RSEM (or salmon/kalisto) for quantifying the abundance of all transcripts in each treatment: I used the Trinity script to call the programs (each has to be installed individually).
- something to visualize genes that is significantly differentially expressed in between treatments (R but not very knowledgeable of this).

## General Steps
1) Clean raw data (Cutadapt, Trimmomatic, Trim_galore)
   * Remove adapters and primers (if present)
   * Remove low-quality sequences
   * Trimmomatic can be run as a part of Trinity!

2) Map reads to either a genome or a transcriptome
   * **No Ref Genome**
     * If no reference is available, assemble a reference transcriptome from all available data: i.e. pool all samples and make a *metatranscriptome*
       * *Trinity* (optionally RNA-spades)
       * Possible things to be aware of:
       * allelic variations, alternative splicing and isoforms
       * Mixed reads if the transcriptome contains reads from hosts and parasite/virus
   *  **With Ref Genome**
      *  Either map to the genomic sequences, or to the predicted genes.
3) Quantify reads
   * Count the number of mapped reads.
     * RPM, RPKM, FPKM, TPM, TMM, DESeq, SCnorm, GeTMM, ComBat-Seq etc....
     * https://www.reneshbedre.com/blog/expression_units.html
     * Normalization depends on what you want to quantify
     * Between samples vs. between genes in the same sample
     * https://github.com/zhangyuqing/ComBat-seq
 

# Part II: 28-30. Sept 2022
## Genome annotation of *Hapolina ericina* UIO028
The host genome is sequenced at JGI. The scaffolds have been separated into "main.genome", "prokaryote" and "redundant".
To use the genome for the transcriptome data, it is necessary to do some annotation. 
1) Masking of low-complexity regions, and repetitive regions. 
   * The common tool to use is repeatmodeller. It takes time to run (~48 hours?), and it might not be finished before the workshop starts... 
   * Most recent version on Saga: RepeatModeler/2.0.2a-foss-2020b
   * See script in *scripts/2022* folder
2) Prepare rna-seq for gene-prediction
   * Align one or more of the samples to the genome with STAR (more is probably better).
     * Current version of STAR on Saga is v2.10.4, I've asked for an update to newest.

3) After masking: run gene prediction with rna-seq data. 
   * Braker2: is probably the best option. Can take RNA data mapped to the genome as "hints" (i.e. where genes are located)
     * https://github.com/Gaius-Augustus/BRAKER
     * https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_to_Braker2.html#gsc.tab=0
     * The setup is complicated with braker2. We have found that installing braker2 from GitHub and exporting the path to the executables is working (other solutions might exist)
     * Augustus: nice if your species is part of the list with gene models (part of braker2)
     * GeneMark ES (part of braker)
     * When predicting the genes from the rna-seq and the 
   * Transdecoder is a simple and fast option if you are in a hurry... (probably not accurate)

4) Functional annotation of predicted genes and identification of reads
   * Eggnog: 
   * Interproscan: proteins as input  (saga: InterProScan/5.47-82.0-GCCcore-9.3.0) 
   * MEME (protein input: motif-scanning only)
   * Diamond (super-fast blast-like algorithm)
   * Kraken2 can be used to identify reads. 
# Mapping RNA-seq to a genome: 
1) STAR is a  splice-aware mapper. Suitable for mapping of eukaryotic reads that might contain splice variants. 
   * Can take gtf-files (i.e. annotations of genes) while mapping. Or features in gtf-files can be used post-maping to extract regions of interest
  -   First, build a database of the genome
  -   Map the reads. 

1b) Kallisto (which we used last time), is a k-mer based pseudo-aligner. Needs a transcriptome/CDS or transcripts as reference. Kan be used against the predicted transcripts produced with braker2. 
    
- FeatureCount in Subreads (also in the package *Rsubread*) to get the gene count matrix from bamfiles produced by STAR: https://sourceforge.net/projects/subread/ . Countmatrix can be imported in DESeq2 or EdgeR.

# Building a transcriptome without a genome. 
Pooling *all* reads might not be the best strategy. It will be very tough on resources to make the assembly. It might not even be feasible. 
- (But I have started the job with spades), but then I cleaned firs with fastp
- Per sample is probably better. 
  1) Trinity is the algorithm that is ususally preferred. But it is relatively slow (40+hours? per sample)
  2) Rna-spades is much faster, clocking in at 25 min(!) pr sample.
   *  Tested on Sample_02-c01vb. 
      * *hard_filtered_transcripts.fasta*: sum = 59 Mbp, n = 71497, ave = 827.12, largest = 49414
      * **transdecoder** on *hard_filtered_transcripts*: 58444 aa

- ggf to gtr


# Further things to explore:
- WGCNA: weighted gene co-expression network analysis
  - https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
- Gene set enrichment analysis:  might only work for humans, or model organisms:
  - DAVID https://david.ncifcrf.gov/