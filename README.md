# TNSeq Analysis Pipeline


## General Steps for Quantifying Barcodes

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