# Title     : DESeq analysis of mBARq data (SKroon)
# Created by: ansintsova
# Created on: 21.10.21
#
#install.packages(c('BiocManager','Hmisc'),
#                 dependencies='Depends',
#                repo = "http://cran.us.r-project.org")
#BiocManager::install(c('DESeq2','apeglm'))
#

library(DESeq2)
args <- commandArgs(trailingOnly=TRUE)
# Sample data file
sdf_file <-args[1]
# Counts file
edf_file <- args[2]
design <- 'tissue'

sdf <- read.csv(sdf_file, row.names=1)
edf <- read.csv(edf_file, row.names=1)
dds <- DESeq2::DESeqDataSetFromMatrix(countData=edf, colData=sdf, design=as.formula(paste('~', design)))
dds$tissue <- relevel(dds$tissue, ref = "inoculum")
dds <- DESeq2::DESeq(dds)
feces_res <- DESeq2::results(dds, contrast=c("tissue","feces","inoculum"))
cc_res <- DESeq2::results(dds, contrast=c("tissue","cc","inoculum"))
write.csv(feces_res,
          '/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/08_21/results/skroon/21-10-2021-skroon_feces_results.csv' )
write.csv(cc_res,
          '/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/08_21/results/skroon/21-10-2021-skroon_cc_results.csv', )
