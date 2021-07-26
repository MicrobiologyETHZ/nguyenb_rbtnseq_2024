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
import statsmodels
import scipy
from datetime import date

from quality_control import generate_DE_dataset, get_fitness_results

def get_barcode_fitness(cnts):
    # Calculate a mean value for all inoculum samples
    inoc_df = cnts[[c for c in cnts.columns if 'd0' in c]]
        
    cnts['inoculum'] = inoc_df.mean(axis=1)
    cnts = cnts.dropna(subset=['inoculum'])  # do I have any NAs?
    cnts = cnts[cnts.inoculum > 0]
    # Caclucalte fitness for each barcode
    fitness = cnts.apply(lambda x: 2 ** x / 2 ** cnts['inoculum']).reset_index()
    return fitness


def get_barcode_fitness_by_library(cnt_df, library, good_samples, outdir, filter_below=0):
    lib_df = cnt_df[cnt_df.library == library].copy()
    lib_df = lib_df[~((lib_df.libcnt.isna()) & (lib_df.phenotype.isna()))]
    library_samples = [s for s in good_samples if s in lib_df.sampleID.unique()]
    sdf, edf, design = generate_DE_dataset(lib_df, library_samples, sample_id='sampleID', filter_below=filter_below)

    _, vst_cnts = get_fitness_results(outdir, dnaid=library.replace("_", "-"), experiment='',
                                        sdf=sdf, edf=edf, design=design)
    vst_df = vst_cnts.rename({'Unnamed: 0': 'barcode'}, axis=1).set_index('barcode')
    
    fit_list = []
    for k,v in sdf.groupby(['dnaid', 'experiment']).indices.items():
        df = vst_df.iloc[:, v]
        fdf = get_barcode_fitness(df)
        fit_list.append(fdf.set_index("barcode"))
    barcode_fitness = pd.concat(fit_list, axis=1).drop(['inoculum'], axis=1)
    return vst_df, barcode_fitness


def get_wt_fitness_by_library(fitness, annotation, phenotype='wt'):
    wt_fitness = fitness.reset_index().rename({'index':'barcode'}, axis=1)
    wt_fitness = wt_fitness.merge(annotation[['barcode', 'phenotype', 'conc']].dropna(), on='barcode')
    wt_fitness = wt_fitness[wt_fitness.phenotype == phenotype]
    return wt_fitness



def get_gene_fitness_by_library(fitness, annotation):
    sample_map = {s:['median'] for s in fitness}
    fitness = fitness.reset_index().rename({'index':'barcode'}, axis=1)  
    fitness = fitness.merge(annotation[['barcode', 'ShortName']], on='barcode')
    gene_fitness = fitness.groupby('ShortName').agg(sample_map).reset_index()
    gene_fitness.columns = [c[0] for c in gene_fitness.columns]
    return gene_fitness


def melt_sampleID(df, idVar=['ShortName'], value_name='fitness'):
    mdf = df.melt(id_vars=idVar, var_name='sampleID', value_name=value_name)
    new = mdf.sampleID.str.split("_", expand=True)
    new.columns = ['mouse', 'day', 'dnaid', 'experiment']
    return pd.concat([mdf, new], axis=1)


def gene_ranksums(gene_values, wt_values):
    return ranksums(gene_values, wt_values)[1]


def fdr_correction(pvals):
    return statsmodels.stats.multitest.multipletests(pvals, alpha=0.05, method='fdr_bh')[1]


def get_library_results(meltGeneFit, meltWtFit, library):
    day_results = []
    for day in meltGeneFit.day.unique():
        if day == 'd0':
            continue
        print(day)
        # Subset by day
        wt_values = meltWtFit[meltWtFit.day == day].groupby('sampleID').fitness.median().values
        day_fitness = meltGeneFit[meltGeneFit.day == day]
        # For each gene run ranksums
        pvals = day_fitness.groupby('ShortName').fitness.apply(gene_ranksums, wt_values=wt_values)
        padj = fdr_correction(pvals.values)
        day_results.append(pd.DataFrame([pvals.values, padj], columns=pvals.index, index=['pval', 'padj']).T.assign(day=day)) 
    results = pd.concat(day_results).assign(library=library)
    results = (meltGeneFit.groupby(['ShortName', 'day']).fitness.median()
               .reset_index().merge(results.reset_index(), on=['ShortName', 'day']))
    return results


def get_library_results_ci(meltGeneCI, ssa_ci, library):
    day_results = []
    for day in meltGeneCI.day.unique():
        if day == 'd0':
            continue
        print(day)
        # Subset by day
        day_ci = meltGeneCI[meltGeneCI.day == day]
        day_ssa = ssa_ci[[i for i in day_ci.sampleID.unique()]]
        # For each gene run ranksums
        pvals = day_ci.groupby('ShortName').ci.apply(gene_ranksums, wt_values=day_ssa)
        padj = fdr_correction(pvals.values)
        day_results.append(pd.DataFrame([pvals.values, padj], columns=pvals.index, index=['pval', 'padj']).T.assign(day=day)) 
    results = pd.concat(day_results).assign(library=library)
    results = (meltGeneCI.groupby(['ShortName', 'day']).ci.median()
               .reset_index().merge(results.reset_index(), on=['ShortName', 'day']))
    return results



# def method2_analysis2(vst, annotation_df, good_samples, sample_id='sampleID',
#                       days=['_d1', '_d2', '_d3', '_d4'], hits=0.05):
#     gene_fitness_dfs = []
#     ci_dfs = []
#     results_dfs = []
#     control_dfs = []
#     barcode_fitness = get_barcode_fitness(vst)
#     for day in days:
#         print(day)
#         controls = get_wt_fitness(barcode_fitness, day, good_samples)
#         gene_fitness = get_gene_fitness(barcode_fitness, day, good_samples)
#         gene_ci, results, controls = get_results(gene_fitness, controls, day, good_samples, hits=hits)

#         if all([d.empty for d in [gene_fitness, gene_ci, results, controls]]):
#             continue
#         df = gene_fitness.reset_index().melt(id_vars=['ShortName'], var_name=sample_id, value_name='Fitness')
#         df['day'] = day.strip('_')
#         gene_fitness_dfs.append(df)
#         df2 = gene_ci.reset_index().melt(id_vars=['ShortName'], var_name=sample_id, value_name='CI')
#         df2['day'] = day.strip('_')
#         ci_dfs.append(df2)
#         results['day'] = day.strip("_")
#         results_dfs.append(results)
#         control_dfs.append(controls)

#     gene_fitness_df = pd.concat(gene_fitness_dfs)
#     ci_df = pd.concat(ci_dfs)
#     results_df = pd.concat(results_dfs)
#     control_df = pd.concat(control_dfs)
#     return barcode_fitness, gene_fitness_df, ci_df, results_df, control_df
