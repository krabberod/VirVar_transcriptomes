# See: http://bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")

library("Rsubread")
fc <- featureCounts(files = "star_bam/Aligned.out.bam", 
                    annot.ext = "star_bam/PkV-RF01_final.gtf", 
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE)
# You will need a list of the gene names, here we are using the same names
# as already used in the DESeq2 anlysis. 
rownames(fc$counts) <- rownames(dds)

# when you have multiple bam in a folder (mock-data): 

files <- list.files("star_bam/", pattern = "bam", full.names =TRUE)
fc <- featureCounts(files,
                    annot.ext = "star_bam/PkV-RF01_final.gtf", 
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE)

rownames(fc$counts) <- rownames(dds)


