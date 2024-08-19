# Salmonella Typhimurium screen identifies shifts in mixed acid fermentation during gut colonization

## Analysis with mBARq

- We have set up small snakemake pipeline to use for data analysis with [mBARq](https://github.com/MicrobiologyETHZ/mbarq).


1. Pre-process Illumina short read data using [standard workflow](https://methods-in-microbiomics.readthedocs.io/en/latest/preprocessing/preprocessing.html).

    
2. Map libraries to reference genome using [mbarq map](https://mbarq.readthedocs.io/en/latest/mapping.html) command. 

```bash
rbseq map -c rbseq_workflow/configs/28-06-23-mapping-config.yaml

```

3. Count the barcodes across samples using [mbarq count](https://mbarq.readthedocs.io/en/latest/counting.html) command.

```
rbseq merge -c rbseq_workflow/configs/03-08-23-counting-config.yaml 

```
    
4. Perform differential abundance analysis using [mbarq analyze](https://mbarq.readthedocs.io/en/latest/analysis.html) command.  

```
rbseq analyze -c rbseq_workflow/configs/03-08-23-counting-config.yaml

```
5. Explore. All of the data generated at each step of the worklow is available for download and further exploration through [mBARq app](https://mbarq.microbiomics.io/).


## Data Analysis

- All of the code used to generate figures for the paper can be found in `code/notebooks/08-23-rbtnseq-analysis.ipynb`.
- The data needed to run the notebook is in `code/notebooks/nguyen.tar.gz` archive. Run `tar -xzvf code/notebooks/nguyen.tar.gz` to extract the archive and make sure the `root` in `code/notebooks/nguyenb_config.yaml` is pointing to the extracted folder.
- Packages required to run the notebook are listed in `code/notebooks/env.yaml`. 
