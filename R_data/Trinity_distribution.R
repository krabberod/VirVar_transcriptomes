# Script for ploting distribution of coverage and length of transcripts in Trinity
# Based on the output of Trinity that 
library(ggplot2)
library(magrittr)

rnaspades_hard_v <- read.table("/Users/anderkkr/Downloads/tmp/Samples_rnaspades_hard_v.fasta.header", header = F)
 
colnames(rnaspades_hard_v) <- c("Sample","node", "length", "cov","gene", "isoform")
str(rnaspades_hard_v)
sort(rnaspades_hard_v$cov, decreasing = TRUE)
hist(rnaspades_hard_v$cov)
# Some with extreme coverage. At least some of them are SSU + LSU
hist(log(rnaspades_hard_v$cov))
# plot(rnaspades_hard_c$Cov, rnaspades_hard_c$length)

rnaspades_hard_v$length

rnaspades_hard_v %>% 
  ggplot(aes(x=log(cov))) + 
  geom_histogram() #+ facet_grid(. ~ Sample)

rnaspades_hard_v %>% 
  ggplot(aes(length)) + 
  geom_histogram() #+ 
  #facet_wrap(. ~ Sample, ncol=5)

unique(rnaspades_hard_v$Sample)

rnaspades_hard_v[rnaspades_hard_v$Sample %in% "Sample_44-c08vb",] %>% 
  ggplot(aes(length)) + 
  geom_histogram() + 
  facet_wrap(. ~ Sample, ncol=5)

rnaspades_hard_v[rnaspades_hard_v$Sample %in% "Sample_44-c08vb",] %>% 
  ggplot(aes(x=log(cov))) + 
  geom_histogram()
