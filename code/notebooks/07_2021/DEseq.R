# Title     : TODO
# Objective : TODO
# Created by: ansintsova
# Created on: 23.01.21
#
#install.packages(c('BiocManager','Hmisc'),
#                 dependencies='Depends',
#                repo = "http://cran.us.r-project.org")
#BiocManager::install(c('DESeq2','apeglm'))
#


library(DESeq2)

args <- commandArgs(trailingOnly=TRUE)

sdf_file <-args[1]

edf_file <- args[2]

design <- args[3]

experiment <- unlist(strsplit(basename(sdf_file), '_sdf'))[[1]]
prefix <- unlist(strsplit(sdf_file, '_sdf'))[[1]]
sdf <- read.csv(sdf_file, row.names=1)
edf <- read.csv(edf_file, row.names=1)


calculate_fitness <- function(sdf, edf, prefix){

    if(design == 'day'){
        print('Design: day')
        dds <- DESeq2::DESeqDataSetFromMatrix(countData=edf, colData=sdf, design=~day)
    }
    else{
        print('Design: experiment + day')
        dds <- DESeq2::DESeqDataSetFromMatrix(countData=edf, colData=sdf, design=~experiment + day)
    }
    dds <- DESeq2::DESeq(dds)
    vsd <- DESeq2::vst(dds, blind=TRUE, nsub=min(1000, nrow(dds)))
    fitness <- list()
    for(day in levels(dds$day)[-1]){
        #fitness[[day]] <- as.data.frame(DESeq2::lfcShrink(dds, coef=paste0("day_", day, "_vs_d0"), type="apeglm"))
        fitness[[day]] <- as.data.frame(DESeq2::results(dds, name=paste0("day_", day, "_vs_d0")))
    }
    write.csv(assay(vsd), paste0(prefix, "_vst.csv"))
    lapply(levels(dds$day)[-1], function(x) write.table(fitness[[x]], paste(prefix, '_fitness_', x, ".txt", sep="")))
}

calculate_fitness(sdf, edf, prefix)

