# Salmonella SL1344 TNSeq Analysis Pipeline

1. Set up small snakemake pipeline to use for data analysis with [tnseq2 package](https://github.com/MicrobiologyETHZ/nccr_tnseq2).

### Mapping of transposon libraries to the SL1344 genomes

1. Pre-process Illumina short read data using standard workflow

    - Create `tnseq-mapping_sample.csv`
    - Preprocess using standard pipeline
    
```bash
cd code/tnseq_snakemake_workflow
# Creates sample sheet for all the libraries to be mapped
python preprocessing/scripts/fastq_dir_to_samplesheet.py -c configs/18-08-21-mapping-config.yaml
# Preprocess (still troubleshooting snakemake pipeline, issues with qout/qerr because of import of workflow)
snakemake --use-conda -k --cluster "qsub -S /bin/bash -V -cwd -o preprocess.qoutfile -e preprocess.q
errfile -pe smp 4 -l h_vmem=8000M" --configfile /nfs/nas22/fs2202/biol_micro_bioinf_nccr
/hardt/nguyenb/tnseq/code/tnseq_snakemake_workflow/configs/18-08-21-mapping-config.yaml
-p -j 1 --max-jobs-per-second 1 preprocess_preprocess

```

2. Map libraries to reference genomes

```bash
snakemake --use-conda -k --cluster "DIR=\$(dirname {params.qoutfile}); mkdir -p \"\${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M" --configfile /nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/code/tnseq_snakemake_workflow/configs/18-08-21-mapping-config.yaml -p -j 10 --max-jobs-per-second 1 map

```



3. (Demultiplex) and Quantify barcodes 


    
 
    
    
    
-  tnseq2 command for mapping:

```
tnseq2 maplib ...
```
- Conda environment used for mapping is in `envs/map.yaml`
- Snakemake command to run this on samples specified in `samples.txt`

```
snakemake --configfile ...
```

### General Steps for Quantifying Barcodes

1. Demultiplex
2. Quantify
3. Run merge_counts.py 
4. Calculate and graph R for all the controls
5. Filter count data based on read distribution in the inoculum
6. Calculate fitness using DESeq2 (for each barcode)
7. Compare fitness of each gene to that of wt using z-scores
8. Visualize the results

Final output table should contain for each gene:


## Results
Two tables:
###Fitness table

Gene ShortName, # of barcodes, mean fitness, sd fitness, num of samples used to calculate fitness 

## Z-score table
z-score, padj, CI for every day. 