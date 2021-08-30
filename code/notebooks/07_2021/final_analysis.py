import collections
import pandas as pd
from pathlib import Path
import datetime as dt
from scipy import stats
from scipy.stats import ranksums
import numpy as np
import math
import subprocess
from scipy.stats import norm
#import statsmodels
from  statsmodels.stats.multitest import multipletests
import statsmodels.stats.multitest

import scipy
from datetime import date


from quality_control import generate_DE_dataset, get_fitness_results


def sigma(lfcSE):
    return np.sqrt(lfcSE.pow(2).sum()) / len(lfcSE)


def calculate_2dist_zscore(u1, s1, u2, s2):
    return (u1 - u2) / np.sqrt((s1 ** 2) + (s2 ** 2))


def to_list(x):
    bc_list = list(x)
    if len(bc_list) == 1:
        return bc_list[0]
    return ", ".join(list(set(x)))


def calculate_comparisons2(fitness, df, control_file):
    """

    fitness: DESeq2 output, log2FoldChange value for each barcode comparing each time point with inoculum
    df: df for 1 experiment and 1 dnaid
    controls: control meta df?
    """
    days = sorted(list(fitness['day'].unique()))
    # days.remove('d0')
    controls = pd.read_table(control_file, names=['barcode', 'phenotype', 'conc'])

    controls = pd.read_table(control_file, names=['barcode', 'phenotype', 'conc'])
    controls['CntrlName']= controls['phenotype'] + controls['conc'].astype(str)
    controls_bc = controls[controls.phenotype == 'wt'].CntrlName.values
    cntrl_df = fitness[fitness.barcode.isin(controls_bc)]
    gene_df = fitness[~fitness.barcode.isin(controls_bc)].rename({'barcode':'ShortName'}, axis=1)
    gene_mean = gene_df.groupby(['ShortName', 'day']).agg(
            {'log2FoldChange': ['mean', 'median'], 'lfcSE': [sigma]}).reset_index()
    gene_mean.columns = ['gene', 'day', 'gene_FC', 'gene_FC_median', 'sigma']
    cntrl_mean = cntrl_df.groupby(['day']).agg({'log2FoldChange': ['mean', 'median'], 'lfcSE': [sigma]})
    cntrl_mean.columns = ['cntrl_FC', 'cntrl_FC_median', 'cntrl_sigma']
    cntrl_mean = cntrl_mean.reset_index()
    gene_mean = gene_mean.merge(cntrl_mean, how='left', on='day')
    gene_mean['zscore'] = gene_mean.apply(
            lambda x: calculate_2dist_zscore(x['gene_FC'], x['sigma'], x['cntrl_FC'], x['cntrl_sigma']), axis=1)
    gene_mean['ci'] = gene_mean.apply(lambda x: 2 ** x['gene_FC'] / 2 ** x['cntrl_FC'], axis=1)
    gene_mean = gene_mean[['gene', 'day', 'gene_FC',  'sigma', 'zscore', 'ci']]
    results = gene_mean.copy()
    results['pval'] = results.zscore.apply(lambda x: scipy.stats.norm.sf(abs(x)) * 2)
    results['padj'] = results.groupby('day').pval.transform(lambda x: multipletests(x, alpha=0.05, method='fdr_bh')[1])
    return results


def get_control_fitness(fitness, control_file):
    controls = pd.read_table(control_file, names=['barcode', 'phenotype', 'conc'])
    cnrls = controls.merge(fitness, how='left', on='barcode').assign(control='yes')
    cnrls = cnrls[['barcode', 'phenotype', 'conc', 'log2FoldChange', 'n_samples', 'day', 'control']]
    cnrls = cnrls[cnrls.day.notnull()]
    days = list(cnrls.day.unique())
    cnrls['gene'] = cnrls['phenotype'] + "-" + cnrls['conc'].astype(str)
    mean_fit = cnrls.groupby(['gene', 'day']).agg({'log2FoldChange':['mean', 'median'], 'barcode':['nunique', to_list]}).reset_index()
    mean_fit.columns = ['gene', 'day', 'gene_FC', 'gene_FC_median', 'num_barcodes', 'barcode']
    mean_fit = mean_fit.merge(cnrls[[ 'gene', 'n_samples', 'day', 'control']], how='left', on=['gene','day'])
    return mean_fit


def final_fitness_table(fitness, exp_df, control_file, results):
    barcode_info = exp_df[['barcode', 'locus_tag', 'ShortName', 'library']].drop_duplicates()
    fitness = fitness.merge(barcode_info, how='left', on=['barcode'])
    bc_per_gene = fitness.groupby(['library', 'ShortName']).agg({'barcode': ['nunique', to_list]}).reset_index()
    bc_per_gene.columns = ['library', 'ShortName', 'num_barcodes', 'barcode']
    num_samples = fitness[['day', 'n_samples']].drop_duplicates()
    control_fit = get_control_fitness(fitness, control_file)

    fit_summary = (results.merge(num_samples, how='left', on='day')
                   .merge(bc_per_gene, how='left', left_on='gene', right_on='ShortName')
                   )
    return pd.concat([fit_summary, control_fit])



def analyze_library2(fdf, sample_id, good_samples, dnaid, experiment, control_file, to_filter, outdir):
    """
     param fdf: already subset datframe
     good_samples: samples passing the correlation cutoff

     If not enough samples, stop the analysis

    # Could add paired option to DESeq2 as now combining results from multiple experiments.

    """

    n_samples = collections.Counter([si.split("_")[1] for si in good_samples])
    d0 = n_samples.pop('d0', 0)
    print(n_samples)
    # Check that for at least one condition have replicates, if not quit
    if (not d0 >= 1) or (not any([i > 1 for i in n_samples.values()])):
        print('oi')
        return pd.DataFrame()
    else:
        # Filter
        print("Filtering Dataset")
        sdf, edf, design = generate_DE_dataset(fdf, good_samples, filter_below=to_filter, sample_id=sample_id)
        # Run DESeq2
        print("Running DESeq2")
        fitness, vst_cnt = get_fitness_results(outdir, dnaid, experiment, sdf, edf, design)
        # Calculate z-scores
        print('Calculating z-scores')
        results = calculate_comparisons2(fitness, fdf, control_file)
        return fitness, results, vst_cnt




