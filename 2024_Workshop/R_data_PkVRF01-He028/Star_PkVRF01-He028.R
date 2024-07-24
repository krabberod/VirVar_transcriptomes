library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggforce)
library(ggrepel)
library("ggbeeswarm")
library(gridExtra)


#### Read Files ####
# Define the path to the directory containing your files
# path_to_files <- "path/to/your/files"
path_to_files <- "ReadsPer/"

# List all the files matching your pattern
files <- list.files(path = path_to_files, pattern = "*.tab", full.names = TRUE)

# Function to read a file and extract the sample name
read_sample_file <- function(file) {
  # Extract sample name include PkVRF01-He028 if you want
  sample_name <- str_extract(basename(file), "Sample_\\d{2}-[a-zA-Z0-9]+")
  # Read the file with appropriate column names
  df <- read_tsv(file, col_names = c("Gene", sample_name), col_types = cols(.default = "c"))
  df[1:2]
  df[[sample_name]] <- as.numeric(df[[sample_name]])
  return(df)
}

# Read all files and store them in a list
list_of_dfs <- lapply(files, read_sample_file)

# Reduce the list of data frames into a single data frame by merging on the "Gene" column
merged_df <- purrr::reduce(list_of_dfs, full_join, by = "Gene")

# Convert the result to a tibble
result_tibble <- as_tibble(merged_df)

# View the resulting tibble
print(result_tibble)

# Find Genes with row sum equal to 0
result_tibble %>%
  filter(rowSums(select(., -Gene)) == 0) # %>% write_csv("Genes_0_count.csv")


# Extract genes that start with Pk (i.e. Virus genes)
Virus <- result_tibble %>%
  filter(str_detect(Gene, "^Pk"))
# Check if any are 0 in read count 
Virus %>%
  filter(rowSums(select(., -Gene)) == 0)

# Extract genes from result_tibble that do not start with Pk or N 
Host <- result_tibble %>%
  filter(!str_detect(Gene, "^Pk")) %>%
  filter(!str_detect(Gene, "^N"))

# add He028_ as prefix to gene name: 
Host$Gene <- paste("He028_", Host$Gene, sep = "")

# Remove genes with 0 read count
# It might be wise to remove gene with low count, not only those that are zero
Host <- Host %>%
  filter(rowSums(select(., -Gene)) > 3)


# Write Host data frame to a CSV file
write_csv(Host, "Host_Genes.csv")
#Write Virus data frame to a CSV file
write_csv(Virus, "Virus_Genes.csv")

# Sort Host by rowsum.
Host <- Host[order(rowSums(select(Host, -Gene)), decreasing = TRUE),]

# Simple plotting of raw data, no normalization: 
Host %>% .[1:16,] %>%
  pivot_longer(cols = -Gene) %>%
  ggplot(aes(x = name, y = value)) +
  geom_col(fill = "skyblue", color = "black") +
  labs(title = "Count of Host Genes in Sample 12",
       x = "Sample",
       y = "Count") +
  theme_minimal() +
  facet_wrap(~Gene, nrow = 4, ncol = 4, scales = "free_y")

# Simple plotting of virus data
# sort virus by rowsum
Virus <- Virus[order(rowSums(select(Virus, -Gene)), decreasing = TRUE),]
green_palette <- brewer.pal(9, "Greens")
Virus %>% .[1:16,] %>%
  pivot_longer(cols = -c(Gene)) %>%
  ggplot(aes(x = name, y = value)) +
  geom_col(fill = green_palette[5], color = "black") +
  labs(title = "Count of Host Genes in Sample 12",
       x = "Sample",
       y = "Count") +
  theme_minimal() +
  facet_wrap(~Gene, nrow = 4, ncol = 4, scales = "free_y")


#### Differential gene expression analysis ####
# Most of the code below is from the DESeq2 vignette, which was updater jan 2024:
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
library(readxl)
library(DESeq2)
library("pheatmap")
library("glmpca")
library("apeglm")
library("genefilter")

# Read metadata file
metadata <- read_excel("Metadata_He028_PkVRF01.xlsx") %>% as.data.frame()

# Change name in the treatment column to infected and control
metadata <- metadata %>%
  mutate(Treatment = ifelse(Treatment == "v", "infected", "control"))

# Extract count data for the Host genes
count_data <- Host %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

all(rownames(metadata) == colnames(count_data))

# Create a DESeq2 dataset
dds_org <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, 
                              design=~Treatment+Replicates+Timepoints)

### Perform DESeq2 analysis ###
dds_org <- DESeq(dds_org)
results <- results(dds_org)

# Explrore the results
nrow(dds_org)

# From the DeSEQ2 vignette:
#> We perform pre-filtering to keep only rows that have a count of at least 10 
#> for a minimal number of samples. The count of 10 is a reasonable choice for 
#> bulk RNA-seq. A recommendation for the minimal number of samples is to 
#> specify the smallest group size, e.g. here there are 4 samples in each group.smallestGroupSize <- 4
smallestGroupSize <- 4
keep <- rowSums(counts(dds_org) >= 10) >= smallestGroupSize
dds_org <- dds_org[keep,]
nrow(dds_org)
rm(keep)

# Variance stabilizing transformation (VST)
vsd_org <- vst(dds_org, blind = FALSE)
dds_org <- estimateSizeFactors(dds_org)

sampleDists <- dist(t(assay(vsd_org)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd_org$dex, vsd_org$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

plotPCA(vsd_org, intgroup = "Timepoints")
plotPCA(vsd_org, intgroup = "Treatment")
plotPCA(vsd_org, intgroup = "Replicates")
plotPCA(vsd_org, intgroup = c("Replicates",  "Timepoints"))

# The ar timepoints are strange... 
### Filter out samples with timepoint "ar" ###
filtered_metadata <- colData(dds_org) %>%
  as.data.frame() %>%
  filter(Timepoints != "ar")

# Subset the count data to include only the filtered samples
filtered_count_data <- count_data[, rownames(filtered_metadata)]

# Create a new DESeq2 dataset with the filtered data
dds <- DESeqDataSetFromMatrix(countData = filtered_count_data, 
                                       colData = filtered_metadata, 
                                       design=~Treatment+Replicates+Timepoints)

# Now proceed with the DESeq2 analysis on the filtered dataset
dds <- DESeq(dds)
# Do the filtering again:
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
rm(keep)

# Variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)
dds <- estimateSizeFactors(dds)
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


plotPCA(vsd, intgroup = "Timepoints")
plotPCA(vsd, intgroup = "Treatment")
plotPCA(vsd, intgroup = "Replicates")

# GLM-NMDS
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Treatment <- dds$Treatment
gpca.dat$Timepoints <- dds$Timepoints
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Timepoints, shape = Treatment)) +
  geom_mark_ellipse(aes(fill = Timepoints), alpha = 0.2) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

# Classic MDS 
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Timepoints, shape = Treatment)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")


pca_data <- plotPCA(vsd, intgroup = "Timepoints", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoints)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
res <- results(dds)
summary(res)
mcols(res, use.names = TRUE)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
# More stringent: 
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Treatment"))

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Treatment", "Timepoints"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Treatment, y = count, color = Timepoints, group = Timepoints)) +
  scale_y_log10() +  geom_beeswarm(cex = 1) + geom_line()

resultsNames(dds)
res <- lfcShrink(dds, coef="Treatment_infected_vs_control", type="apeglm")
plotMA(res, ylim = c(-5, 5))

# Gene clustering
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

anno <- as.data.frame(colData(vsd)[, c("Treatment", "Timepoints")])
pheatmap(mat, annotation_col = anno, cluster_cols =F)

# Subset data for 'infected' group
infected_data <- vsd[, vsd$Treatment == "infected"]
topVarGenes_infected <- head(order(rowVars(assay(infected_data)), decreasing = TRUE), 100)
mat_infected <- assay(infected_data)[topVarGenes_infected, ]
mat_infected <- mat_infected - rowMeans(mat_infected)
anno_infected <- as.data.frame(colData(infected_data)[, c("Treatment", "Timepoints")])
heatmap_infected <- pheatmap(mat_infected, annotation_col = anno_infected, cluster_cols = T,cluster_rows = T, plot = FALSE)

# Subset data for 'control' group
control_data <- vsd[, vsd$Treatment == "control"]
topVarGenes_control <- head(order(rowVars(assay(control_data)), decreasing = TRUE), 100)
mat_control <- assay(control_data)[topVarGenes_control, ]
mat_control <- mat_control - rowMeans(mat_control)
anno_control <- as.data.frame(colData(control_data)[, c("Treatment", "Timepoints")])
heatmap_control <- pheatmap(mat_control, annotation_col = anno_control, cluster_cols = FALSE, cluster_rows = T, plot = FALSE)

# Plot the heatmaps together in different panels
grid.arrange(heatmap_infected$gtable, heatmap_control$gtable, ncol = 2)


##### Differential expression in time series ####
# Time series example: 
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments

# Define the reference level for Timepoints
dds_filtered$Timepoints <- relevel(dds_filtered$Timepoints, ref = "0")


dds_interaction <- DESeqDataSet(dds, ~ Timepoints + Treatment + Timepoints:Treatment)
print(table(metadata$Timepoints, metadata$Treatment))


# Trying to include the interaction term in the design formula
# However timepoint 17 lacks the infected stage and 34 the control stage
# They are removed from the anlaysis

# Filter out samples with Timepoints "ar", "17", and "34"
filtered_metadata <- colData(dds) %>%
  as.data.frame() %>%
  filter(Timepoints != "ar" & Timepoints != 17 & Timepoints != 34)

filtered_metadata$Timepoints <- factor(filtered_metadata$Timepoints)
filtered_metadata$Treatment <- factor(filtered_metadata$Treatment)

# Drop unused factor levels
filtered_metadata$Timepoints <- droplevels(filtered_metadata$Timepoints)
filtered_metadata$Treatment <- droplevels(filtered_metadata$Treatment)

# Check the levels of Timepoints and Treatment
print(levels(filtered_metadata$Timepoints))
print(levels(filtered_metadata$Treatment))
filtered_count_data <- count_data[, rownames(filtered_metadata)]
print(table(filtered_metadata$Timepoints, filtered_metadata$Treatment))

dds_interaction <- DESeqDataSetFromMatrix(countData = filtered_count_data, 
                                          colData = filtered_metadata, 
                                          design = ~ Timepoints + Treatment + Timepoints:Treatment)

# From the Tutorial: 
#> The following chunk of code performs a likelihood ratio test, 
#> where we remove the treatment-specific differences over time. 
#> Genes with small p values from this test are those which at one or more time 
#> points after time 0 showed a treatment-specific effect. Note therefore that 
#> this will not give small p values to genes that moved up or down over time 
#> in the same way in both strains.
smallestGroupSize <- 4
keep <- rowSums(counts(dds_interaction) >= 10) >= smallestGroupSize
dds_interaction <- dds_interaction[keep,]
nrow(dds_interaction)
rm(keep)

dds_interaction <- DESeq(dds_interaction, test="LRT", 
                         reduced = ~ Treatment + Timepoints)
dds_interaction <- nbinomLRT(dds_interaction,reduced = ~ Treatment + Timepoints, maxit=3000)

res_interaction <- results(dds_interaction)
res_interaction$symbol <- mcols(res_interaction)$symbol


head(res_interaction[order(res_interaction$padj),], 10)
infection <- plotCounts(dds_interaction, which.min(res_interaction$padj), 
                   intgroup = c("Timepoints","Treatment"), returnData = TRUE)
infection$Timepoints <- as.numeric(as.character(infection$Timepoints))

ggplot(infection,
       aes(x = Timepoints, y = count, color = Treatment)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()


# Plot with wesanderson colors
  ggplot(infection,
         aes(x = Timepoints, y = count, color = Treatment)) + 
    geom_point() + 
    stat_summary(fun.y=mean, geom="line") +
    scale_y_log10() +
    scale_color_manual(values = wes_palette("AsteroidCity2"))

resultsNames(dds_interaction)
res02 <- results(dds_interaction, name="Timepoints02.Treatmentinfected", test="Wald")
res02[which.min(res02$padj),]

res10 <- results(dds_interaction, name="Timepoints10.Treatmentinfected", test="Wald")
res10[which.min(res10$padj),]

# We can furthermore cluster significant genes by their profiles. We extract a 
# matrix of the log2 fold changes using the coef function. Note that these are 
# the maximum likelihood estimates (MLE). For shrunken LFC, one must obtain 
# them one coefficient at a time using lfcShrink.


betas <- coef(dds_interaction)
colnames(betas)
topGenes <- head(order(res_interaction$padj),50)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

# Timepoints_10_vs_0 etc are baseline
# Timepoints02.Treatmentinfected is the interaction term, the reaction to infection!

# Make a summary of the results
summary(res_interaction)


