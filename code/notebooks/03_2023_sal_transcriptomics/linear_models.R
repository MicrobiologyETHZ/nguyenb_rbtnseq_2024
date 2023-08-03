# Installation of the package for https://academic.oup.com/bioinformatics/article/37/Supplement_1/i34/6319701#409224236
# In terminal loaded R module and opened R
#library(devtools)
#install_github('biobakery/MTX_model')

library(tidyverse)
data_dir <- "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/03_23_transcriptomics/modeling/synth_mgx_mtx"

input_data <- read.csv(file.path(data_dir, "true-exp.mtx_abunds_count.csv"), row.names=1)
meta_data <- read.csv(file.path(data_dir, "true-exp.mtx_abunds_meta.csv"), row.names=1)
covariate_data <- read.csv(file.path(data_dir, "true-exp.bug_abunds_edited.csv"), row.names=1)

library(MTXmodel)

fit_data <- MTXmodel(t(input_data), meta_data, 'test_out', 
transform = 'LOG', 
fixed_effects = c('Phenotype'),
reference = "Phenotype,0",
standardize = FALSE,
input_dnadata = t(covariate_data),
rna_dna_flt="strict")



fit_data <- MTXmodel(t(input_data), meta_data, 'test_out_rna', 
transform = 'LOG', 
fixed_effects = c('Phenotype'),
standardize = FALSE,
)

###
df <- read.csv("/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/03_23_transcriptomics/modeling/test_df.csv")

m <- glm("tss_log ~ Phenotype", data =  df)
lm_summary <- summary(m)$coefficients
para <- as.data.frame(lm_summary)[-1, -3]
para$name <- rownames(lm_summary)[-1]
para
