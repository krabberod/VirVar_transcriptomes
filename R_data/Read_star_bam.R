# See: http://bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html
# On saga: 
# Ask for interactive computing node:
#
# module load R/4.1.2-foss-2021b


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Rsubread")

library("Rsubread")
  fc <- featureCounts(files = "star_bam/Sample1.bam", 
                    annot.ext = "star_bam/PkV-RF01_final.gtf", 
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE)
# You will need a list of the gene names, here we are using the same names
# as already used in the DESeq2 analysis. 
rownames(fc$counts) <- rownames(dds)

# when you have multiple bam in a folder (mock-data): 

# files <- list.files("star_bam/Sample_01-a01va/", pattern = "bam", full.names =TRUE)
files<-list.dirs("star_bam/") %>% .[-1] %>% file.path(.,"Aligned.out.bam")

fc <- featureCounts(files,
                    annot.ext = "star_bam/PkV-RF01_final.gtf", 
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE)

rownames(fc$counts) <- rownames(dds)

colnames(fc$counts) <- samples[1:9]

dds<-DESeqDataSetFromMatrix(fc$counts, 
                              colData = metadata[1:9,], 
                              design=~Treatment+Replicates+Timepoints)

