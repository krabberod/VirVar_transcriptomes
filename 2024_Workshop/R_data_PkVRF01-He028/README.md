# Scripts for deferential Gene Expression Analysis in R. 

This repository contains the scripts for deferential Gene Expression Analysis in R, based on the data from the transcriptome analysis of the PkVRF01-He028 interaction.

The input data is the tabulated ReadsPerGene.out.tab files from Star alignment of RNA-seq to the reference genome of PkVRF01 and He028. The data is processed using the R package DESeq2.

The metadata file is a tabulated file with the following columns (in excel format):
| Sample_name    | Sample | Timepoints | Treatment | Replicates |
|----------------|--------|------------|-----------|------------|
| Sample_01-a01va | a01va  | 0          | v         | a          |
| Sample_02-a01vb | a01vb  | 0          | v         | b          |
| Sample_03-a01vc | a01vc  | 0          | v         | c          |
| Sample_04-a01ca | a01ca  | 0          | c         | a          |
| Sample_05-a01cb | a01cb  | 0          | c         | b          |
| Sample_06-a01cc | a01cc  | 0          | c         | c          |
| Sample_07-a02vb | a02vb  | 02         | v         | b          |
| Sample_08-a02vc | a02vc  | 02         | v         | c          |
| Sample_09-a02ca | a02ca  | 02         | c         | a          |
| Sample_10-a02cb | a02cb  | 02         | c         | b          |

The main R script is [Star_PkVRF01-He028.R](./Star_PkVRF01-He028.R).

# Data Analysis Workflow

## Libraries
The analysis utilizes several R libraries for data manipulation, visualization, and statistical analysis, including `tidyverse`, `ggplot2`, `RColorBrewer`, `ggforce`, `ggrepel`, `ggbeeswarm`, `gridExtra`, `readxl`, `DESeq2`, `pheatmap`, `glmpca`, and `apeglm`.

## Reading and Merging Data
1. **File Reading**: Files with a specific pattern (`*.tab`) are read from a designated directory.
2. **Extracting Sample Names**: Sample names are extracted from file names.
3. **Data Import**: Each file is read into a data frame with appropriate column names and types.
4. **Merging Data**: All data frames are merged into a single data frame based on the "Gene" column.

## Data Cleaning and Filtering
1. **Identify Zero Sum Genes**: Genes with a total read count of zero across all samples are identified.
2. **Classifying Genes**:
   - **Virus Genes**: Genes starting with "Pk" are classified as virus genes.
   - **Host Genes**: Genes not starting with "Pk" or "N" are classified as host genes.
3. **Prefix Addition**: "He028_" is prefixed to host gene names.
4. **Filter Low Count Genes**: Host genes with low read counts (≤3) are removed.

## Data Export
- Host and virus gene data frames are exported to CSV files.

## Data Visualization
1. **Host Genes**: Raw counts of the top 16 host genes are plotted using `ggplot2`.
2. **Virus Genes**: Raw counts of the top 16 virus genes are plotted using a green palette.

## Differential Gene Expression Analysis
1. **Metadata Import**: Metadata is read from an Excel file.
2. **DESeq2 Dataset Creation**: A DESeq2 dataset is created using the host gene count data and metadata.
3. **DESeq2 Analysis**: Differential expression analysis is performed.
4. **Filtering**: Low count genes are filtered out, and variance stabilizing transformation (VST) is applied.
5. **Distance Matrix and Clustering**: Sample distances are calculated and visualized using heatmaps and PCA plots.
6. **Filtering Specific Timepoints**: Samples with specific timepoints are filtered out, and a new DESeq2 dataset is created and analyzed.

## Advanced Analysis
1. **GLM-NMDS and MDS**: Generalized PCA and classic MDS are performed and visualized.
2. **PCA**: Principal Component Analysis (PCA) is performed and visualized for different groups.
3. **Differential Expression Results**: Results are summarized, and significant genes are identified and visualized.
4. **Gene Clustering**: Significant genes are clustered based on their profiles and visualized using heatmaps.

## Time Series Analysis
1. **Interaction Analysis**: Interaction terms for timepoints and treatment are included in the DESeq2 design.
2. **Likelihood Ratio Test (LRT)**: An LRT is performed to identify genes with treatment-specific effects over time.
3. **Clustering**: Significant genes are clustered by their profiles using log2 fold changes, and heatmaps are generated.

## Summary of Results
The results of the differential expression and time series analyses are summarized, highlighting the key findings and significant genes.
