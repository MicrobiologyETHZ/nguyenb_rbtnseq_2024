library(dada2)
library(tidyverse)


asvs <- read.table("/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENERAL_METAB_ANALYSIS_PAN/data/processed/IMB/NGUY23-1/NGUY23-1.asvs.tsv", sep='\t', header=TRUE)
seqs <- asvs$seq
rownames(asvs) <- asvs$seq
asvs <- asvs %>% select(-tax, -seq)
asvs <- t(asvs)
head(asvs[,1:2])

taxa <- assignTaxonomy(seqs, "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/03_23_transcriptomics/ref/silva_nr99_v138_train_set.fa", multithread=TRUE)
