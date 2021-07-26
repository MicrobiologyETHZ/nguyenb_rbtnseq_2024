import pandas as pd
from pathlib import Path
import sys


def read_tnseq_count_file(f):
    df = pd.read_table(f, sep='\t', header=None)
    df.columns = ['barcode', 'cnt', 'position', 'seq', 'strand']
    df['codeID'] = f.stem.split(".")[0].split("_")[1]
    return df


def merge_counts(count_dir, meta_dir, dnaid):
    # concat all the count files
    counts = [read_tnseq_count_file(f) for f in (Path(count_dir) / dnaid).iterdir()]
    df = pd.concat(counts)
    df['dnaid'] = dnaid

    # map columns using the meta file
    meta_file = (pd.read_table(Path(meta_dir) / f"{dnaid}_metadata.txt", header=None,
                               names=['codeID', 'library', 'experiment', 'DN1', 'mouse', 'day', 'organ'],
                               dtype={'codeID': str})
                 .drop('DN1', axis=1))

    df = df.merge(meta_file, on='codeID')
    df['sampleID'] = df['mouse'] + "_" + df['day']

    return df


def get_library(mapping_file):
    library = pd.read_table(mapping_file, header=None, index_col=0, names=['lib_counts', 'position', 'chr',
                                                                          'strand', 'norm_count', 'locus', 'gene'])
    library['library'] = Path(mapping_file).parent.stem
    return library

def all_libraries(map_dir):
    files = [c/"barcode_map.txt" for c in Path(map_dir).iterdir()]
    df_list = [get_library(f) for f in files]
    return pd.concat(df_list)


def add_mapping_info(merged_df, mapping_df):
    mdf = mapping_df.reset_index().rename({'index':'barcode'}, axis=1)[['barcode', 'locus', 'gene', 'library']]
    fdf = merged_df.merge(mdf, how = 'left', left_on =['barcode', 'library'], right_on = ['barcode', 'library'])
    fdf.to_csv(Path(count_dir) / f"{dnaid}/{dnaid}_merged_counts.csv")
    return fdf



if __name__ == "__main__":
    count_dir = sys.argv[1]
    meta_dir = sys.argv[2]
    dnaid = sys.argv[3]
    map_dir = sys.argv[4]

    df = merge_counts(count_dir, meta_dir, dnaid)
    fdf = add_mapping_info(df, all_libraries(map_dir))
    print(fdf.head())