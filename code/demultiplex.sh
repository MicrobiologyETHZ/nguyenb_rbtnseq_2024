#$ -cwd
#$ -S /bin/bash
#$ -N tnseq_demux
#$ -V
#$ -pe smp 32
#$ -l h_vmem=4G
#$ -e logs/demux.error.log
#$ -o logs/demux.out.log

ml Python

tnseq demux -bc /nfs/home/ansintsova/multiplex_codes.txt -r1 /nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/scratch/chris/hardt/nguyenb/tnseq_sample_2023/data/dnaid2023_BKDL202599908-1A_HCWCNCCX2_L2_1.fq.gz -s 0 -o /science/ansintsova/nguyenb/demux