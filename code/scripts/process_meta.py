"""
Assume metadata comes in the following format:

Sample
1315-10-library11_1-TV3379-inoculum
1457-6-library11_2-TV3522A-w437-d1-feces
sequencing run-sample barcode-mutant library#-mouse experiment number
(which can be subdivided e.g. TV5585A, TV5585B, TV5585C, as I can use different libraries
for the same batch of animals--in other words treat subdivisions as independent experiments )
-animal number (inoculum or unenriched_inoculum)-time point (d1, d2, d3, d4)-content

"""
import pandas as pd


def read_metadata(file):
    df = pd.read_table(file)
    samples = df.columns
    df = df.replace(regex={"inoculum$": 'inoculum-d0-inoculum', "unenriched$": "unenriched_inoculum-d0-inoculum"})

    for sample in samples:
        print(sample)
        sdf = df[sample].dropna()
        print(sdf.head())
        sdf = sdf.str.split('-', expand=True)
        sdf.columns = ['seqRun', 'demuxCode', 'library', 'batch', 'mouse', 'day', 'tissue']
        sdf.to_csv(f'../data/metadata/{sample}_metadata.txt', sep='\t',
                   columns=['demuxCode', 'library', 'batch', 'batch', 'mouse', 'day', 'tissue'],
                   header=None, index=None)
    return None

if __name__ == "__main__":
    file = "../../data/Untreated-LCM-Sequencing-libraries.txt"
    read_metadata(file)