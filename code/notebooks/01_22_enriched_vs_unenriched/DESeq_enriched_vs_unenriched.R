# Title     : DESeq analysis of mBARq data (B. Nguyen)
# Created by: ansintsova
# Created on: 05.01.22
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
library <- args[3]
design <- args[4] # options are "mouse" or "mouse + experiment"
outDir <- args[5]

sdf <- read.csv(sdf_file, row.names=1)
edf <- read.csv(edf_file, row.names=1)
dds <- DESeq2::DESeqDataSetFromMatrix(countData=edf, colData=sdf, design=as.formula(paste('~', design)))
dds$mouse <- relevel(dds$mouse, ref = "unenriched_inoculum")
dds <- DESeq2::DESeq(dds)
results = DESeq2::results(dds, contrast=c("mouse", "inoculum", "unenriched_inoculum"))
write.table(results, paste0(outDir, '/', library, '-results-inoculum', ".csv"))
