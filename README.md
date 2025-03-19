Amazon molly loss of heterozygosity analysis using whole genome resequencing data from Warren et al. (2018). 

Dataset includes 19 hybrid Poecilia formosa samples, 5 parental P. latipinna samples, 4 parental P. mexicana samples, and outgroup samples from P. velifera and P. gillii. Pipeline is implemented in snakemake form. 

To recreate, please obtain the reference genome at (INSERT RERFERENCE SOURCE HERE) and modify the data/refInfo.tsv file with the appropriate path info.

This pipeline is designed to be run in a computing environment with SLURM job submission. The snakemake command should be run in a se