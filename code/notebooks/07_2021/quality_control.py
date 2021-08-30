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


def load_files(samples, results_dir):
    df_list = []
    cntrl_list = []
    for sample in samples:
        df_list.append(pd.read_csv(Path(results_dir)/sample/"merged_counts.csv", index_col=0))
        cntrl_list.append(pd.read_csv(Path(results_dir)/sample/"merged_controls.csv", index_col=0))
    return pd.concat([pd.concat(df_list), pd.concat(cntrl_list)])


def calculate_correlation(exp_df, control_file, for_each='sampleID', how='log', cutoff=0.8):
    """
    Subset counts for control barcodes
    Calculate correlation on log counts (log), log counts, but keep 0 (log_w_0), or raw data (raw)

    """
    controls = pd.read_table(control_file, names=['barcode', 'phenotype', 'conc'])
    control_cnts = exp_df[exp_df.barcode.isin(controls.barcode)].copy()
    if how == 'raw':
        col1 = 'conc'
        col2 = 'cnt'
    else:
        control_cnts['logConc'] = np.log10(control_cnts['conc'])
        col1 = 'logConc'
        if how == 'log':
            control_cnts['logCnts'] = np.log10(control_cnts['cnt'])
        elif how == 'log_w_0':
            control_cnts['logCnts'] = np.log10(control_cnts['cnt'].replace({0: 1}))
        col2 = 'logCnts'
    corr_df = control_cnts.groupby(['phenotype', for_each])[[col1, col2]].corr()
    corr_df = corr_df.reset_index()
    corr_df = corr_df[corr_df['level_2'] == col1].drop(['level_2', col1], axis=1)
    corr_df.columns = ['phenotype', 'sampleID', 'R']

    good_samples = corr_df[(corr_df.R**2 > cutoff) & (corr_df.phenotype == 'wt')].sampleID.values
    return corr_df, good_samples


def filter_inoculum(exp_df, filter_below=0, sample_id='sampleID'):
    filt_df = exp_df.copy()
    if 'ShortName' in filt_df.columns:
        filt_df = filt_df.drop(['ShortName'], axis=1)
    if 'locus_tag' in filt_df.columns:
        filt_df = filt_df.drop(['locus_tag'], axis=1)
    filt_df = (filt_df.drop_duplicates()
               .pivot(index='barcode', columns=sample_id, values='cnt'))

    filt_df = filt_df.fillna(0)
    columns_to_filter = [f for f in filt_df.columns if 'inoculum' in f]

    filt_df = filt_df[(filt_df[columns_to_filter] >= filter_below).all(1)]
    return filt_df


def filter_samples(exp_df, good_samples, sample_id='sampleID'):
    return exp_df.copy()[exp_df[sample_id].isin(good_samples)]


def run_command(args):
    """Run command, transfer stdout/stderr"""
    result = subprocess.run(args)
    try:
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        raise e
        
        
def generate_DE_dataset(exp_df, good_samples, sample_id='sampleID', filter_below=0):

    sample_data = exp_df[[sample_id, 'mouse', 'day', 'tissue', 'dnaid', 'experiment']].set_index(sample_id).drop_duplicates()
    sample_data = sample_data.loc[sample_data.index.intersection(good_samples)]
    sample_data, design = check_sdf(sample_data)
    expr_data = filter_samples(exp_df, good_samples, sample_id)

    expr_data = filter_inoculum(expr_data, filter_below=filter_below, sample_id=sample_id)

    expr_data = expr_data[list(sample_data.index)].reset_index()
    print(f"Number of unique experiments after filtering: {sample_data.experiment.nunique()}")
    print(f"Design: {design}")
    return sample_data, expr_data, design


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