import collections
import pandas as pd
from pathlib import Path
from scipy import stats
import numpy as np
import math
import subprocess
from scipy.stats import norm
import statsmodels
import scipy
from datetime import date
import os


def load_samples(mapped_files,  metadata_df, unmapped_files, control_file=""):
    """
    Take a list of samples with mapped barcodes, unmapped_barcodes, return 1 dataframe for mapped barcodes
    and a control barcode df if control_file is provided
    :param mapped_files:
    :param unmapped_files:
    :param control_file:
    :return:
    """
    df = pd.concat([pd.read_csv(f, index_col=0).assign(sampleName=f.name.split("_counts")[0]) for f in mapped_files])
    df = df.merge(metadata_df, on='sampleName', how='left')
    if control_file:
        control_df = pd.read_table(control_file, header=None, index_col=0)
        control_df.columns = ['barcode', 'phenotype', 'conc']
        df_unmapped = pd.concat([(pd.read_csv(f, index_col=0)
                                  .assign(sampleName=f.name.split("_counts")[0])) for f in unmapped_files])
        control_df = control_df.merge(df_unmapped, on='barcode', how='left')
        control_df = control_df.merge(metadata_df, on='sampleName', how='left')
    else:
        control_df = pd.DataFrame()
    return df, control_df



def get_control_df():
    pass


def calculate_correlation(control_df, concentration_col='conc',
                          cnt_col='cnt', phenotype_col='phenotype',
                          for_each='sampleID', how='log', wt_phenotype='wt', cutoff=0.8):
    """
    Given a control_df:

            barcode	    phenotype	conc	  barcode_cnt	    sampleID
0	TACCCAGAGCACACTCA	wt	        0.0015	  10034.0	        dnaid2030_3
1	TACCCAGAGCACACTCA	wt	        0.0015	  6340.0	        dnaid2030_6
2	TACCCAGAGCACACTCA	wt	        0.0015	  7938.0	        dnaid2030_4

    Calculate correlation on log counts (log), log counts, but keep 0 (log_w_0), or raw data (raw)
    Has to have a column with concentrations for each barcode and with count for each barcode,
    and a column phenotype for each barcode (ex. 'wt')
    """
    if how == 'raw':
        col1 = concentration_col
        col2 = cnt_col
    else:
        control_df['logConc'] = np.log10(control_df[concentration_col])
        col1 = 'logConc'
        if how == 'log':
            control_df['logCnts'] = np.log10(control_df[cnt_col])
        elif how == 'log_w_0':
            control_df['logCnts'] = np.log10(control_df[cnt_col].replace({0: 1}))
        col2 = 'logCnts'
    num_observations = control_df.groupby([for_each, phenotype_col])[col2].count().reset_index()
    num_observations.columns = [for_each, phenotype_col, 'num_observations']
    corr_df = control_df.groupby([phenotype_col, for_each])[[col1, col2]].corr()
    corr_df = corr_df.reset_index()
    corr_df = corr_df[corr_df['level_2'] == col1].drop(['level_2', col1], axis=1)
    corr_df.columns = [phenotype_col, for_each, 'R']
    corr_df = corr_df.merge(num_observations, on=[for_each, phenotype_col])
    fdf = corr_df[corr_df.num_observations > 9]
    good_samples = fdf[(fdf.R**2 > cutoff) & (fdf.phenotype == wt_phenotype)][for_each].values
    return corr_df, good_samples


# def filter_inoculum(exp_df, filter_below=0, sample_id='sampleID'):
#     filt_df = exp_df.copy()
#     if 'ShortName' in filt_df.columns:
#         filt_df = filt_df.drop(['ShortName'], axis=1)
#     if 'locus_tag' in filt_df.columns:
#         filt_df = filt_df.drop(['locus_tag'], axis=1)
#     filt_df = (filt_df.drop_duplicates()
#                .pivot(index='barcode', columns=sample_id, values='cnt'))
#
#     filt_df = filt_df.fillna(0)
#     columns_to_filter = [f for f in filt_df.columns if 'inoculum' in f]
#
#     filt_df = filt_df[(filt_df[columns_to_filter] >= filter_below).all(1)]
#     return filt_df


# def filter_samples(exp_df, good_samples, sample_id='sampleID'):
#     return exp_df.copy()[exp_df[sample_id].isin(good_samples)]


def run_command(args):
    """Run command, transfer stdout/stderr"""
    result = subprocess.run(args)
    try:
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        raise e
        
        
def generate_DE_dataset(exp_df, good_samples, sample_id='sampleID', filter_below=0):
    pass

    # sample_data = exp_df[[sample_id, 'mouse', 'day', 'tissue', 'dnaid', 'experiment']].set_index(sample_id).drop_duplicates()
    # sample_data = sample_data.loc[sample_data.index.intersection(good_samples)]
    # sample_data, design = check_sdf(sample_data)
    # expr_data = filter_samples(exp_df, good_samples, sample_id)
    #
    # expr_data = filter_inoculum(expr_data, filter_below=filter_below, sample_id=sample_id)
    #
    # expr_data = expr_data[list(sample_data.index)].reset_index()
    # print(f"Number of unique experiments after filtering: {sample_data.experiment.nunique()}")
    # print(f"Design: {design}")
    # return sample_data, expr_data, design


def check_sdf(sdf):
    

    # Check more than 1 sample per day, if not drop day
    mice_per_day = sdf[sdf.day != 'd0'].groupby('day').mouse.nunique()
    good_days = list(mice_per_day[mice_per_day > 1].index) + ['d0']
    sdf = sdf[sdf.day.isin(good_days)]
    
    # Check number of experiments for design
    if sdf.experiment.nunique() == 1:
        design = "day"
    else:
        design = "experiment+day"
    return sdf, design


def get_fitness_results(fitness_dir, dnaid, experiment, sdf, edf, design):
    sdf_path = Path(fitness_dir) / f"{dnaid}_{experiment}_sdf.csv"
    edf_path = Path(fitness_dir) / f"{dnaid}_{experiment}_edf.csv"
    sdf.to_csv(sdf_path)
    edf.set_index('barcode').to_csv(edf_path)
    #rpath = Path(__file__).parent.absolute()
    #rpath = Path('/Users/ansintsova/git_repos/tnseq2.0/tnseq2/src')
    r = run_command(['Rscript', './DEseq.R', sdf_path, edf_path, design])
    fitness = pd.concat(
        [pd.read_table(f, sep=' ').assign(day=f.stem.split("_")[3]) for f in Path(fitness_dir).iterdir() if
         f"{dnaid}_{experiment}_fitness" in f.stem])
    vst_counts = pd.read_csv(Path(fitness_dir)/f"{dnaid}_{experiment}_vst.csv")
    for file in [f for f in Path(fitness_dir).iterdir() if f"{dnaid}_{experiment}_fitness" in f.stem]:
        os.remove(file)
    os.remove(sdf_path)
    os.remove(edf_path)
    os.remove(Path(fitness_dir)/f"{dnaid}_{experiment}_vst.csv")
    
    n_samples = sdf.groupby('day').mouse.nunique().to_dict()
    fitness['n_samples'] = fitness.day.map(n_samples)
    fitness = fitness.reset_index().rename({'index': 'barcode'}, axis=1)
    return fitness, vst_counts


def get_skewed(wt_fitness):
    skewed = (abs(np.log2(wt_fitness.drop(['conc'], axis=1).median())) > 1)
    skewed = skewed[skewed].index.values
    return skewed