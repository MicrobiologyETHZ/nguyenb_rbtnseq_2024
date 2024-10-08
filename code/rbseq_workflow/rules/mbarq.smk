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
    sampleData = pd.read_csv(config['sampleData'])
    library = sampleData.loc[sampleData.sampleID == wildcards.sample, 'library'].values[0]
    return Path(config['mappingDir'])/f'{library}/{library}.annotated.csv'


rule quantify_one:
    input: DATADIR / 'demux/{sample}.fasta.gz'
    output: OUTDIR / 'counts/{sample}_mbarq_counts.csv'
    params:
        barcode_map=get_mapping_file,
        outdir=lambda wildcards: OUTDIR /'counts',
        transposon= config['tn'],
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.quant.qout',
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.quant.qerr',
        prefix=lambda wildcards: f'{wildcards.sample}',
        scratch=500,
        mem=8000,
        time=235
    log:
        log=OUTDIR / 'logs/{sample}.quant.log'
    conda:
        'mbarq_lite'
    threads:
        8
    shell: 'mbarq count -f {input} -m {params.barcode_map} '
           ' -o {params.outdir} -n {params.prefix} -tn {params.transposon} &> {log.log}'


def get_lib_samples(wildcards):
    sampleData = pd.read_csv(config['sampleData'])
    sampleData = sampleData[sampleData.dnaid != 'dnaid1315']
    samples = sampleData.loc[sampleData.library == wildcards.lib, 'sampleID'].unique()
    return [OUTDIR / f'counts/{sample}_mbarq_counts.csv' for sample in samples]

def get_lib_samples_as_string(wildcards):
    sampleData = pd.read_csv(config['sampleData'])
    sampleData = sampleData[sampleData.dnaid != 'dnaid1315']
    samples = sampleData.loc[sampleData.library == wildcards.lib, 'sampleID'].unique()
    return ",".join([str(OUTDIR / f'counts/{sample}_mbarq_counts.csv') for sample in samples])

rule merge:
    input: get_lib_samples
    output: OUTDIR/'counts/{lib}_mbarq_merged_counts.csv'
    params:
        count_files = get_lib_samples_as_string,
        out_dir = OUTDIR/'counts',
        qoutfile = lambda wildcards: OUTDIR /f'logs/{wildcards.lib}.merge_counts.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/{wildcards.lib}.merge_counts.qerr',
        name = '{lib}',
        scratch = 500,
        mem = 8000,
        time = 235
    conda:
        'mbarq_lite'
    log:
        log = OUTDIR /'logs/{lib}.merge_counts.log'
    threads:
        8
    shell:
        "mbarq merge -i {params.count_files} -o {params.out_dir} "
        "-a Name  -n {params.name} &> {log.log}"


rule analyze_one:
    input:
        OUTDIR/'counts/{lib}_mbarq_merged_counts.csv'
    output:
        OUTDIR/'analysis/{lib}_rra_results.csv'
    params:
        sample_data=lambda wildcards:Path(config["metaDir"])/f"{wildcards.lib}_sample_data.csv",
        out_dir=OUTDIR / 'analysis',
        control_file = config['controlStrains'],
        baseline = config['baseline'],
        treatment_col = config['treatment_col'],
        norm_method = config['normMethod'],
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.lib}.analysis.qout',
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.lib}.analysis.qerr',
        name='{lib}',
        scratch=500,
        mem=8000,
        time=235
    conda:
        'mbarq_lite'
    log:
        log = OUTDIR /'logs/{lib}.analysis.log'
    threads:
        8
    shell:
        "mbarq analyze -i {input} -o {params.out_dir} -s {params.sample_data} "
        "-c {params.control_file} --baseline {params.baseline} -n {params.name} "
        "--treatment_column {params.treatment_col} -g locus_tag --norm_method {params.norm_method} --filter_low_counts 100"
        " &> {log.log}"


rule analyze_inoc:
    input:
        OUTDIR/'counts/{lib}_mbarq_merged_counts.csv'
    output:
        OUTDIR/'analysis/{lib}_unenriched_inoculum_rra_results.csv'
    params:
        sample_data=lambda wildcards:Path(config["metaDir"])/f"{wildcards.lib}_sample_data.csv",
        out_dir=OUTDIR / 'analysis',
        control_file = config['controlStrains'],
        norm_method = config['normMethod'],
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.lib}_inoc.analysis.qout',
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.lib}_inoc.analysis.qerr',
        name= lambda wildcards: f'{wildcards.lib}_unenriched_inoculum',
        scratch=500,
        mem=8000,
        time=235
    conda:
        'mbarq_lite'
    log:
        log = OUTDIR /'logs/{lib}_inoc.analysis.log'
    threads:
        8
    shell:
        "mbarq analyze -i {input} -o {params.out_dir} -s {params.sample_data} "
        "-c {params.control_file} --baseline unenriched_inoculum -n {params.name} "
        "--treatment_column mouse  -g Name --norm_method {params.norm_method} --filter_low_counts "
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
