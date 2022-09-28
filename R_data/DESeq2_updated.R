### Firs load packages. You might have to install some of them. 
# 
# Exampels: 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

# These are the libraries used: 
library("tximport")
library("DESeq2")
library("tximportData")
library("tidyverse")
library("ggplot2")
library("gplots")
library("vsn")
library("IHW")
library("readxl")
library("rhdf5")
library("apeglm")

# dir <- system.file("extdata", package = "tximportData")
# samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
# files <- file.path(dir, "kallisto_boot", samples$run, "abundance.h5")

# Import Kallisto-results 
files<-list.dirs("kallisko_results") %>% .[-1] %>% file.path(.,"abundance.h5")

samples <- str_remove(files, "/abundance.h5") %>% str_remove("kallisko_results/")
names(files) <-samples
txi.kallisto <- tximport(files, type = "kallisto", 
                         txOut = TRUE, ignoreAfterBar = T)
metadata<-readxl::read_xlsx("Metadata.xlsx")

# Quick heatmap to explore the data:
# based on counts and , wihtout reordering the samples:
heatmap(txi.kallisto$counts)
heatmap(txi.kallisto$abundance, Colv=NA,  Rowv = NA)


# Create the DESeq object: 
dds<-DESeqDataSetFromTximport(txi.kallisto, 
                              colData = metadata, 
                              design=~Treatment+Replicates+Timepoints)
dds <- DESeq(dds)                                 

# Examples from the tutorial
resultsNames(dds)
res <- results(dds, name = "Timepoints_13_vs_01")
resLFC <- lfcShrink(dds, coef="Timepoints_13_vs_01", type="apeglm")
resLFC

# We can order our results table by the smallest p value:
resOrdered <- resLFC[order(resLFC$pvalue),]

# How many adjusted p-values were less than 0.01?
sum(resLFC$padj < 0.01, na.rm=TRUE)

# See what can be done with the results() function: 
?results
  
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)


# Independent hypothesis weighting
# (unevaluated code chunk)

resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.05, na.rm=TRUE)
metadata(resIHW)$ihwResult

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# Plotting the gene with lowest p-value:  
plotCounts(dds, gene=which.min(res$padj), 
           intgroup="Timepoints")

p1 <- plotCounts(dds, 522, intgroup=c("Timepoints","Treatment","Replicates"),
                 returnData=TRUE)


d <- plotCounts(dds, gene=which.min(res$padj), 
                intgroup="Timepoints", 
                returnData=TRUE)


# Added treatment as color
library("ggplot2")
library("ggeasy")
library("ggpubr")
ggplot(p1 , aes(x=Timepoints, y=count, color = factor(Treatment))) + 
  geom_point(position=position_jitter(w=0.1,h=0.0)) +
  easy_add_legend_title("Treatment") +
  ggtitle(paste("Gene"))
  # + 
  #scale_y_log10(breaks=c(25,100,400))


ordered_genes <- rownames(resOrdered) %>% 
  str_split_fixed("_",2) %>% .[,2] %>%  as.numeric()
genes_of_interest <- ordered_genes[1:9]

list_of_plots <-list()

for(i in 1:length(genes_of_interest))
  {
  list_of_plots[[i]]  <- plotCounts(dds, genes_of_interest[i], 
                intgroup=c("Timepoints","Treatment","Replicates"),
                 returnData=TRUE) %>% 
    ggplot(aes(x=Timepoints, y=count, color = factor(Treatment))) + 
    geom_point(position=position_jitter(w=0.1,h=0.0)) +
    easy_add_legend_title("Treatment") +
    ggtitle(paste("Gene",  genes_of_interest[i]))
  }

ggarrange(plotlist=list_of_plots, ncol = 3, nrow = 3)

colData(dds)
resMFType <- results(dds, contrast=c("Timepoints","02","ar"))


# Transformation: 
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Heatmap with clustering and aggregation of rows: 
# The function also allows to aggregate the rows using kmeans clustering. 
# This is advisable if number of rows is so big that R cannot handle their 
# hierarchical clustering anymore, roughly more than 1000. 

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)[,c("Replicates","Timepoints")])

pheatmap(assay(ntd)[select,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


plotPCA(ntd, intgroup="Timepoints",ntop = 500)






### Redo without the control. I.e. Subset to only virus

dds_vir <- dds[, dds$Treatment %in% "va"]
# If you want to remove the re-infected
`%notin%` <- Negate(`%in%`)
dds_vir<-dds_vir[, dds_vir$Timepoints %notin% "ar"]

resultsNames(dds_vir)
res_vir <- results(dds_vir, name = "Timepoints_13_vs_01")
resLFC_vir <- lfcShrink(dds_vir, coef="Timepoints_13_vs_01", type="apeglm")
resLFC_vir

plotCounts(dds_vir, gene=which.min(res_vir$padj), intgroup="Timepoints")
d <- plotCounts(dds_vir, gene=which.min(res_vir$padj), intgroup="Timepoints", 
                returnData=TRUE)

ggplot(d, aes(x=Timepoints, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0.0)) + 
  scale_y_log10(breaks=c(25,100,400))


ntd_vir <- normTransform(dds_vir)
meanSdPlot(assay(ntd_vir))
pheatmap(assay(ntd_vir)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

plotPCA(ntd_vir, intgroup="Timepoints",ntop = 500)


#etc. etc 