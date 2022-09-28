# Script for ploting distribution of coverage and length of transcripts in Trinity
# Based on the output of Trinity that 
library(ggplot2)
library(magrittr)

rnaspades_hard_c <- read.table("trinity/Samples_rnaspades_hard_c.fasta.header", header = F)
 
colnames(rnaspades_hard_c) <- c("Sample","node", "length", "cov","gene", "isoform")
str(rnaspades_hard_c)
sort(rnaspades_hard_c$cov, decreasing = TRUE)
hist(rnaspades_hard_c$cov)
# Some with extreme coverage. At least some of them are SSU + LSU
hist(log(rnaspades_hard_c$cov))
# plot(rnaspades_hard_c$Cov, rnaspades_hard_c$length)

rnaspades_hard_c$length

rnaspades_hard_c %>% 
  ggplot(aes(x=log(cov))) + 
  geom_histogram() #+ facet_grid(. ~ Sample)

rnaspades_hard_c %>% 
  ggplot(aes(length)) + 
  geom_histogram() + 
  facet_wrap(. ~ Sample, ncol=5)

