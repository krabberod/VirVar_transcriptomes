library(tidyverse)
library(ggplot2)

# Script for plotting gene expression of multiple genes in the same plot
# Takes a DeSeq2 object (dds), extracts the normalised counts
# Rarranges data and plots 
# Needs a list of "genes of interest", if you have not made one yet 
# you can make one using the number of the gene:
# genes_of_interest <- c(38,637,892,940)

# Here I use the genes_of_interest from the examples script "DESeq2_updated.R":

# Make sure that some gene have been selected: 
genes_of_interest 


# Extract the Normalized count data for the genes from the dds object
list_of_genes_to_plot <-data.frame(Sample = colnames(dds))
for(i in 1:length(genes_of_interest)){
  tmp <- data.frame()
  tmp <- plotCounts(dds, genes_of_interest[i], 
                  intgroup=c("Timepoints","Treatment","Replicates"),
                  returnData=TRUE, normalized = TRUE ) 
  colnames(tmp) <- c(paste0("gene_",genes_of_interest[i]),"Timepoints", "Treatment", "Replicates")
  tmp$Sample <- colnames(dds)
  list_of_genes_to_plot <- full_join(list_of_genes_to_plot, tmp)
  }


# Rearrange the data to make plotting easier: 
list_of_genes_to_plot_wide <- list_of_genes_to_plot %>% 
  gather(gene, count, grep("gene", colnames(list_of_genes_to_plot)), factor_key=TRUE) 

# Plot. Comment out the second line if want to add the control group:
list_of_genes_to_plot_wide %>%
  filter(Treatment != "c") %>% # Removes the controlgroup
  ggplot(aes(x = Sample, y = count)) +
  geom_line(aes(group = gene), color="grey") + #Connect the dots
  geom_point(aes(color = factor(gene))) +
  easy_add_legend_title("Genes") +
  easy_rotate_x_labels(angle = 90)




