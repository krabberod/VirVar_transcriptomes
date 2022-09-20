# VirVar transcriptomes
Resources for VirVar. Workshop 27-29 Sept 2021


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
   * Trimmomatic can be run as a part of the Trinity!

2) Map reads to either a genome or a transcriptome
   * **No Ref Genome**
     * If no reference is available, assemble a reference transcriptome from all available data: i.e. pool all samples and make a *metatranscpritome*
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
 