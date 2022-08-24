import pandas as pd
from pathlib import Path
# Map

rule sample_mapping:
        input:
            fq1 = DATADIR/"clean_reads/{sample}/{sample}.1.fq.gz",
        output:
            map = OUTDIR/'maps/{sample}/{sample}.annotated.csv',
        params:
            sample = '{sample}',
            outdir = lambda wildcards: OUTDIR/f'maps/{wildcards.sample}',
            transposon = config['tn'],
            genome = config['genome'],
            gff = config['gff'],
            qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.maps.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.maps.qout',
            scratch = 6000,
            mem = 8800,
            time = 1400,
            l = config["map_filter"]
        log:
            log = OUTDIR/'logs/{sample}.map.log',
        conda:
            'mbarq_lite'
        threads:
            16
        shell:
            "mbarq map -f {input.fq1}  -o {params.outdir} -tn {params.transposon} "
            "-l {params.l} -g {params.genome} -n {params.sample} -a {params.gff}  &> {log.log} "


# Count
#
def get_mapping_file(wildcards):
    metadata_file = Path(config['metaDir'])/f'{wildcards.sample}_metadata.edited.txt'
    df = pd.read_table(metadata_file, usecols=[0,1], names=['code', 'library'], dtype={'code':str, 'library':str})
    wc = wildcards.code#.split('_')[1]
    library = df.loc[df.code == wc, 'library'].values[0]
    return Path(config['mappingDir'])/f'{library}/{library}.annotated.csv'


rule quantify_one:
    input: DATADIR / 'demux/{sample}/{sample}_{code}.fasta.gz'
    output: OUTDIR / 'counts/{sample}/{sample}_{code}_mbarq_counts.csv'
    params:
        barcode_map=get_mapping_file,
        outdir=lambda wildcards: OUTDIR / f'counts/{wildcards.sample}',
        transposon= config['tn'],
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}_{wildcards.code}.quant.qout',
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}_{wildcards.code}.quant.qerr',
        prefix=lambda wildcards: f'{wildcards.sample}_{wildcards.code}',
        scratch=500,
        mem=8000,
        time=235
    log:
        log=OUTDIR / 'logs/{sample}_{code}.quant.log'
    conda:
        'mbarq_lite'
    threads:
        8
    shell: 'mbarq count -f {input} -m {params.barcode_map} '
           ' -o {params.outdir} -n {params.prefix} -tn {params.transposon} &> {log.log}'


rule merge:
    input: OUTDIR/'counts/{sample}.done'
    output: OUTDIR/'counts/{sample}_mbarq_merged_counts.csv'
    params:
        count_dir = lambda wildcards: OUTDIR/f'counts/{wildcards.sample}',
        out_dir = OUTDIR/'counts',
        qoutfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.merge_counts.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/{wildcards.sample}.merge_counts.qerr',
        name = '{sample}',
        scratch = 500,
        mem = 8000,
        time = 235
    conda:
        'mbarq_lite'
    log:
        log = OUTDIR /'logs/{sample}.merge_counts.log'
    threads:
        8
    shell:
        "mbarq merge -d {params.count_dir} -o {params.out_dir} "
        "-a Name  -n {params.name} &> {log.log}"


rule analyze:
    input:
        OUTDIR/'counts/{sample}_mbarq_merged_counts.csv'
    output:
        OUTDIR/'analysis/{sample}_rra_results.csv'
    params:
        sample_data=config['sampleData'],
        out_dir=OUTDIR / 'analysis',
        control_file = config['controlStrains'],
        baseline=config['baseline'],
        treatment_col=config['treatmentColumn'],
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.analysis.qout',
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.analysis.qerr',
        name='{sample}',
        scratch=500,
        mem=8000,
        time=235
    conda:
        'mbarq_lite'
    log:
        log = OUTDIR /'logs/{sample}.analysis.log'
    threads:
        8
    shell:
        "mbarq analyze -i {params.input} -o {params.out_dir} "
        "--baseline {params.baseline} "
        "--treatment_column {params.treatment_col} -g Name "
        " &> {log.log}"


rule quantify_db:
    input: DATADIR /'barseq_samples/{sample}.fastq.gz'
    output: DATADIR/ 'counts/{sample}_mbarq_counts.csv'
    params:
        barcode_map=DATADIR/'TnSeq_SB2B_ML5_l5.annotated.csv',
        outdir= DATADIR/'counts',
        transposon=config['tn'],
        qoutfile=lambda wildcards: DATADIR / f'logs/{wildcards.sample}.quant.qout',
        qerrfile=lambda wildcards: DATADIR / f'logs/{wildcards.sample}.quant.qerr',
        scratch=500,
        mem=8000,
        time=235
    log:
        log=DATADIR / 'logs/{sample}.quant.log'
    conda:
        'mbarq_lite'
    threads:
        8
    shell:
        'mbarq count -f {input} -m {params.barcode_map} '
        ' -o {params.outdir}  -tn {params.transposon} --rev_complement -e 0 &> {log.log}'