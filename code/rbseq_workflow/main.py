import argparse
import subprocess
import shlex
import shutil
import os
from pathlib import Path
import click
import yaml

@click.group()
def main():
    pass

from .scripts import fastq_dir_to_samplesheet as fds


# RBSeq

@main.command(help="Generate samplesheet from a directory of FastQ files.")
@click.option('--configfile', '-c', default='', help='Configuration File')
@click.option("-i", "--fastq_dir", help="Folder containing raw FastQ files.")
@click.option("-o", "--sample_file", help="Output samplesheet file.")
@click.option("-r1", "--read1_extension", type=str,  default="_R1.fq.gz",
              help="File extension for read 1.")
@click.option("-r2", "--read2_extension", type=str, default="_R2.fq.gz",
              help="File extension for read 2.")
@click.option("-sn", "--sanitise_name", is_flag=True,
              help="Whether to further sanitise FastQ file name to get sample id. Used in conjunction with "
                   "--sanitise_name_delimiter and --sanitise_name_index.")
@click.option("-sd", "--sanitise_name_delimiter", type=str, default="_",
              help="Delimiter to use to sanitise sample name.", )
@click.option("-si", "--sanitise_name_index", type=int, default=1,
              help="After splitting FastQ file name by --sanitise_name_delimiter "
                   "all elements before this index (1-based) will be joined to create final sample name.",)
def samples(configfile, fastq_dir, sample_file, read2_extension, read1_extension, sanitise_name,
           sanitise_name_delimiter, sanitise_name_index):
    click.echo("Running RBSeq Pipeline")
    if configfile:
        click.echo(f"Config file: {configfile}")
        with open(configfile) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
        fds.fastq_dir_to_samplesheet(
            fastq_dir=config['dataDir'],
            samplesheet_file=config["samples"],
            read1_extension=config['fq_fwd'],
            read2_extension=config['fq_rvr'],
            sanitise_name=config['sanitise_name'],
            sanitise_name_delimiter=config['name_delimiter'],
            sanitise_name_index=config['name_index'],
        )
    else:
        fds.fastq_dir_to_samplesheet(
            fastq_dir=fastq_dir,
            samplesheet_file=sample_file,
            read1_extension=read1_extension,
            read2_extension=read2_extension,
            sanitise_name=sanitise_name,
            sanitise_name_delimiter=sanitise_name_delimiter,
            sanitise_name_index=sanitise_name_index,
        )


@main.command()
@click.option('--config', '-c',  help='Configuration File')
@click.option('--local', is_flag=True, help="Run on local machine")
@click.option('--dry', is_flag=True, help="Show commands without running them")
def clean(config, local, dry):
    click.echo("Running RBSeq Analysis Pipeline: preprocessing sequencing data")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format('locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'preprocess', smk_file, dry, local)
    click.echo(" ".join(cmd))


@main.command()
@click.option('--config', '-c',  help='Configuration File')
@click.option('--local', is_flag=True, help="Run on local machine")
@click.option('--dry', is_flag=True, help="Show commands without running them")
def map(config, local, dry):
    click.echo("Running RBSeq Analysis Pipeline: mapping libraries data")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format('locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'map', smk_file, dry, local)
    click.echo(" ".join(cmd))


@main.command()
@click.option('--config', '-c',  help='Configuration File')
@click.option('--local', is_flag=True, help="Run on local machine")
@click.option('--dry', is_flag=True, help="Show commands without running them")
def merge(config, local, dry):
    click.echo("Running RBSeq Analysis Pipeline: counting and merging samples data")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format('locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'merge_all', smk_file, dry, local)
    click.echo(" ".join(cmd))


def snakemake_cmd(config, analysis, smk_file, dry, local, no_conda=False):
    if dry:
        cmd = shlex.split(f'snakemake -s {smk_file} --configfile {config} -np {analysis} ')
    elif local:
        cmd = shlex.split(f'snakemake -s {smk_file} --configfile {config} --use-conda -j 1 {analysis} ')
    else:
        rstring = r'"DIR=$(dirname {params.qoutfile}); mkdir -p \"${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M"'
        if no_conda:
            part1 = shlex.split(f'snakemake --configfile {config} -s {smk_file} -k --cluster ')
        else:
            part1 = shlex.split(f'snakemake --configfile {config} -s {smk_file} --use-conda -k --cluster ')
        part2 = shlex.split(f'{rstring}')
        part3 = shlex.split(f' -p -j 8 --max-jobs-per-second 1 {analysis}')
        cmd = part1 + part2 + part3
    wdPath = Path(__file__).parent.absolute()
    subprocess.check_call(cmd, cwd=wdPath)
    return cmd


if __name__ == "__main__":
    main()
