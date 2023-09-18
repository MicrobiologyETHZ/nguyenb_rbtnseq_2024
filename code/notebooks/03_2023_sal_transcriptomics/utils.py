

import pandas as pd
from pathlib import Path
import sys
import plotly.express as px
import pyranges as pr
import yaml
import numpy as np

class GenomeAnnot:
    genome_map = {'CP015399.2': 'YL32',
                  'CP015400.2': 'KB18',
                  'CP015401.2': 'I48',
                  'CP015402.2': 'YL27',
                  'CP015403.2': 'YL45',
                  'CP015404.2': 'I46',
                  'CP015405.2': 'YL58',
                  'CP015406.2': 'YL31',
                  'CP015407.2': 'YL2',
                  'CP015408.2': 'I49',
                  'CP015409.2': 'YL44',
                  'CP015410.2': 'KB1',
                  'GCF_000364265': 'ASF519',
                  'FQ312003.1': 'SL1344',
                  'FQ312003.1;FQ312003.1': 'SL1344',
                  'HE654725.1': 'SL1344',
                  'HE654726.1': 'SL1344',
                  'HE654724.1': 'SL1344',
                  'contig_15': 'contig_15',
                  'contig_21': 'contig_21',
                  'contig_26': 'contig_26',
                  'contig_46': 'contig_46',
                  'AQFU02000001.1': 'ASF 502',
                  'AQFU02000002.1': 'ASF 502',
                  'AQFU02000003.1': 'ASF 502',
                  'CP097573.1': 'ASF500',
                  'NZ_CP097810.1': 'ASF356',
                  'NZ_AQFR02000001.1': 'ASF360',
                  'NZ_AQFR02000002.1': 'ASF360',
                  'NZ_AQFR02000003.1': 'ASF360',
                  'NZ_CP097561.1': 'ASF361',
                  'NZ_AQFT02000001.1': 'ASF492',
                  'NZ_AQFT02000002.1': 'ASF492',
                  'NZ_AQFT02000003.1': 'ASF492',
                  'NZ_AQFT02000004.1': 'ASF492',
                  'NZ_AQFT02000005.1': 'ASF492',
                  'NZ_AQFT02000006.1': 'ASF492',
                  'NZ_AQFV02000001.1': 'ASF519',
                  'NZ_AQFV02000002.1': 'ASF519',
                  'NZ_AQFV02000003.1': 'ASF519',
                  'NZ_AQFV02000004.1': 'ASF519',
                  'NZ_AQFV02000005.1': 'ASF519',
                  'NZ_AQFV02000006.1': 'ASF519',
                  'NZ_CP097562.1': 'ASF457'
                  }
    annotation_columns = ['Chromosome',  'Feature', 'Start', 'End', 'Strand', 'ID',
                          'Name', 'locus_tag', 'gene_biotype', 'product']

    def __init__(self, gff_file, ):
        self.gff_file = gff_file
        self.feature = "gene"
        self.annot = self.process_gff()

    def process_gff(self):
        gff = pr.read_gff3(self.gff_file).as_df()[self.annotation_columns]
        return gff[gff['Feature'] == self.feature]

    def annotate_df(self, df):
        fdf = df.merge(self.annot, on='ID', how='left')
        fdf['genome'] = fdf['Chromosome'].replace(
            self.genome_map)
        return fdf


class CountDataSet:
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.count_data = pd.DataFrame()

    def load_count_files(self):
        pass


# Currently not looking at these


class HtseqCounts(CountDataSet):
    count_col = 'count'

    def load_count_files(self):
        files = list(self.data_dir.rglob('*.txt'))
        df_list = []
        for f in files:
            df = pd.read_table(f, names=['Name', 'count'], header=None).assign(
                sample_id=f.stem.split(".")[0])
            df['Name'] = df.Name.str.split("gene-", expand=True)[1]
            df = df.dropna(subset=['Name'])
            df['genome'] = [self.genome_map.get(name.split(
                "_")[0], 'SL1344') for name in df.Name.values]
            df_list.append(df)
        self.count_data = (pd.concat(df_list).rename({self.count_col: 'read_counts',
                                                     'Name': 'locus_tag'}, axis=1)
                           .merge(self.annot, on='locus_tag', how='left'))


class SalmonCounts(CountDataSet):
    count_col = 'NumReads'
    gene_col = 'Name'

    def load_count_files(self):
        files = list(self.data_dir.rglob('quant.sf'))

        df_list = []
        for f in files:
            name = f.parent.stem.split("_quant")[0]
            print(name)
            df = pd.read_table(f).assign(sample_id=name)
            df = df.rename(
                {self.count_col: 'salmon_read_counts'}, axis=1)
            df['ID'] = (df[self.gene_col].str.split('ID=', expand=True)[1]
                        .str.split(";", expand=True)[0])
            df = df.drop(columns=[self.gene_col])
            df_list.append(df)
        self.count_data = pd.concat(df_list)


class FeatureCounts(CountDataSet):
    count_col = None

    def load_count_files(self):
        files = list(self.data_dir.rglob("*.count.txt"))
        df_list = []
        for f in files:
            name = f.stem.split(".count")[0]
            print(name)
            df = pd.read_table(f, comment='#').assign(sample_id=name)
            df.columns = ['ID', 'chr', 'start', 'end',
                          'strand', 'length', 'fc_read_counts', 'sample_id']
            df = df[['ID', 'fc_read_counts', 'sample_id']]
            df_list.append(df)
        self.count_data = pd.concat(df_list)

    @property
    def summary_df(self):
        files = list(self.data_dir.rglob("*.count.txt.summary"))
        df_list = []
        for f in files:
            df = pd.read_table(f)
            name = df.columns[1].split("/")[-1].split('.')[0]
            df = df.assign(sample_id=name)
            df.columns = ['status', 'read_counts', 'sample_id']
            df_list.append(df)
        fdf = pd.concat(df_list)
        summary = fdf.groupby('sample_id').read_counts.sum().reset_index()
        summary.columns = ['sample_id', 'total']
        summary = (summary.merge(fdf[fdf.status == 'Assigned'][['read_counts', 'sample_id']], on='sample_id')
                   .rename({'read_counts': 'assigned'}, axis=1)
                   .merge(fdf[fdf.status == 'Unassigned_Unmapped'][['read_counts', 'sample_id']], on='sample_id')
                   .rename({'read_counts': 'unmapped'}, axis=1)
                   .merge(fdf[fdf.status == 'Unassigned_NoFeatures'][['read_counts', 'sample_id']], on='sample_id')
                   .rename({'read_counts': 'no_feature'}, axis=1))
        summary['percent_assigned'] = summary['assigned']/summary['total']*100
        summary['percent_unmapped'] = summary['unmapped']/summary['total']*100
        summary['percent_no_feature'] = summary['no_feature'] / \
            summary['total']*100
        return summary


class SushiCounts(CountDataSet):
    count_col = "total_insertcount"
    gene_col = "#reference"

    def load_count_files(self):
        files = list(self.data_dir.rglob("*ushicounts"))
        df_list = []
        for f in files:
            name = f.stem.split(".")[0]
            print(name)
            df = pd.read_table(
                f, usecols=[0, 2, 6, 7, 8]).assign(sample_id=name)
            df = df.rename(columns={self.count_col: "sushi_insertcount"})
            df['ID'] = df['#reference'].str.split(
                ';', expand=True)[0].str.split('ID=', expand=True)[1]
            df = df.drop(columns=[self.gene_col])
            df_list.append(df)
        self.count_data = pd.concat(df_list)

        # self.count_data = fdf.merge(self.annot, on='ID', how='left')
        # self.count_data["genome"] = self.count_data['Chromosome'].replace(
        #     self.genome_map)

