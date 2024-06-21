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

