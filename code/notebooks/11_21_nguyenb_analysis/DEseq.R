# Title     : DESeq analysis of mBARq data (B. Nguyen)
# Created by: ansintsova
# Created on: 09.11.21
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
#sdf_file <- "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/08_21/results/nguyenb/nguyen_test_sdf.csv"
# Counts file
edf_file <- args[2]
#edf_file <- "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/08_21/results/nguyenb/nguyen_test_edf.csv"
library <- args[3]
design <- args[4] # options are "day" or "day + experiment"
outDir <- args[5]

sdf <- read.csv(sdf_file, row.names=1)
edf <- read.csv(edf_file, row.names=1)
dds <- DESeq2::DESeqDataSetFromMatrix(countData=edf, colData=sdf, design=as.formula(paste('~', design)))
dds$day <- relevel(dds$day, ref = "d0")
dds <- DESeq2::DESeq(dds)

fitness <- list()
for(day in levels(dds$day)[-1]){
        fitness[[day]] = as.data.frame(DESeq2::results(dds, contrast=c("day", day, 'd0')))
    }
#outdir <- "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/08_21/results/nguyenb/"
lapply(levels(dds$day)[-1], function(x) write.table(fitness[[x]], paste(outDir, '/', library, '-results-', x, ".csv", sep="")))



