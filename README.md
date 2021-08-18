# Salmonella SL1344 TNSeq Analysis Pipeline

1. Set up small snakemake pipeline to use for data analysis with [tnseq2 package](https://github.com/MicrobiologyETHZ/nccr_tnseq2).

### Mapping of transposon libraries to the SL1344 genomes

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