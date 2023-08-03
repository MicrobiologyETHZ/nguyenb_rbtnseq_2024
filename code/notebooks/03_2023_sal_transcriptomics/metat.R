library(DESeq2)
library(edgeR)
library(tidyverse)
source("03_2023_sal_transcriptomics/deseq_fx.R")
library(config)
library(data.table)

# Testing synthetic data
configs <- config::get(file="03_2023_sal_transcriptomics/config.yaml")
root <- configs$root
synth_dir <- file.path(root, configs$synth_dir)
synth_sample_data_file <- file.path(synth_dir, "true-exp-meta.csv")
synth_count_file <- file.path(synth_dir, 'true-exp-counts.csv')
conditions <- "Phenotype"
out_dir <- file.path(synth_dir)
lfct <- 0.0
alpha <- 0.01
date <- Sys.Date()
prefix <- paste0(date, "_true-exp-deseq-")

# First running DESeq2 on the full count data without normalisation. Count data had +1 added to it. Otherwise DESeq2 fails (due to low coverage)
deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix, TRUE, 10)

dds <- deseq_output$dds
sd <- deseq_output$sample_data
pattern1 <- '1'
pattern2 <- '0'
get_all_results(sd, pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)


# Now running within taxon normalisation

count_data <- read.csv(synth_count_file)
count_data <- count_data %>% mutate(genome = str_split(ID, "_", simplify = T)[, 1])
sample_data <- read.csv(synth_sample_data_file)
conditions <- "Phenotype"
prefix <- paste0(date, "_true-exp-deseq-within-taxon-")
deseq_output <- deseq_on_metat_taxon(count_data, sample_data, conditions, "ID", "sample_id", "genome")
dds <- deseq_output$dds
get_all_results(sd, pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)


# Now running on each individual bug one by one. 
prefix <- paste0(date, "_true-exp-deseq-taxon-one-by-one-")
bugs <- count_data$genome %>% unique()
df_list <- list()
for (i in seq_along(bugs)){
    data <- count_data %>% filter(genome == bugs[[i]])
    print(bugs[[i]])
    raw_data <- align_samples_and_counts(data, sample_data, conditions, "ID", "sample_id")
    dds <- DESeqDataSetFromMatrix(countData = raw_data$count_data, colData = raw_data$sample_data, design = ~group)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    dds <- DESeq(dds, quiet = FALSE)
    df_list[[i]] <- results(dds, contrast = c("group", pattern1, pattern2), alpha = alpha, lfcThreshold = lfct) %>% as.data.frame() %>% mutate(bug=bugs[[i]])
    
}

final_result <- rbindlist(df_list, fill=TRUE)
write.csv(final_result, file.path(out_dir, paste0(prefix, pattern1, "_vs_", pattern2, "_l", lfct, "a", alpha, "_results.csv")))

# Now applying this to OLIGO and LCM data
out_dir <- file.path(root, configs$output_dir)

count_data_file <- file.path(root, configs$deseq_count_file)
sample_data_file <- file.path(root, configs$sample_data_file)

count_data <- read.csv(count_data_file)
sample_data <- read.csv(sample_data_file)
conditions <- "Mouse,Treatment"
prefix <- paste0(date, "_oligo-lcm-within-taxon-")
deseq_output <- deseq_on_metat_taxon(count_data, sample_data, conditions, "ID", "sample_id", "genome")
dds <- deseq_output$dds


pattern1 <- "LCM_D"
pattern2 <- "LCM_PBS_D1"
get_all_results(colData(dds), pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)


pattern1 <- "Oligo_PBS"
pattern2 <- "LCM_PBS_D1"
get_all_results(colData(dds), pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)


oligo <- sample_data %>% filter(Mouse == 'Oligo')

oligo_counts <- count_data %>% select(c(c('ID', 'genome'), oligo$sample_id))
oligo_counts <- oligo_counts %>% filter(genome != "ASF457")
conditions <- 'Treatment'
prefix <- paste0(date, "_oligo-alone-within-taxon-")

oligo_output <- deseq_on_metat_taxon(oligo_counts, oligo, conditions, "ID", "sample_id", "genome")
dds <- oligo_output$dds
pattern1 <- "LPS"
pattern2 <- "PBS"
get_all_results(colData(dds), pattern1, pattern2, dds, out_dir, lfct, alpha, prefix)

write.csv(oligo_output$norm_counts, file.path(out_dir, paste0(prefix, 'norm_cnts.csv')))



# have to set the following env variables: export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
# Settings set to active terminal
configs <- config::get(file="03_2023_sal_transcriptomics/config.yaml")
root <- configs$root
sample_data_file <- file.path(root, configs$sample_data_file)

# OLIGO EXPERIMENT
conditions <- "Treatment"
out_dir <- file.path(root, configs$output_dir)

lfct <- 0.0
alpha <- 0.01
date <- Sys.Date()

# Analysing strains together
prefix <- paste0(date, "_oligo-metat-fc-deseq-")
count_file <- file.path(root, configs$oligo_fc_raw)
deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix, TRUE)

dds <- deseq_output$dds
sd <- deseq_output$sample_data

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)

# Analysing different strains seperately
# I48

prefix <- paste0(date, "_oligo-i48-fc-deseq-")
count_file <- file.path(root, configs$oligo_i48_fc_raw)
deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sample_data

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)

# YL27

prefix <- paste0(date, "_oligo-yl27-fc-deseq-")
count_file <- file.path(root, configs$oligo_yl27_fc_raw)

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sample_data

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)


#-----------------------------------------------------------------------


# YL32

prefix <- "08-06-23-YL32-high-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL32_high_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)



# YL32 - all

prefix <- "08-06-23-YL32-all-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL32_all_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)




# YL27

prefix <- "08-06-23-YL27-all-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL27_all_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)



# YL58

prefix <- "08-06-23-YL58-med-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL58_med_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)



# YL58

prefix <- "08-06-23-YL58-all-"
count_file <- "data/05-23_skroon-rnaseq/02-06-23_YL58_all_samples.csv"

deseq_output <- run_deseq_on_featcnts(sample_data_file, count_file, conditions,
                                      out_dir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData

pattern1 <- 'LPS'
pattern2 <- 'PBS'
get_all_results(sd, pattern1, pattern2, dds, out_dir, prefix, annotations)

