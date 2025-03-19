Amazon molly loss of heterozygosity analysis using whole genome resequencing data from Warren et al. (2018). 

Dataset includes 19 hybrid Poecilia formosa samples, 5 parental P. latipinna samples, 4 parental P. mexicana samples, and outgroup samples from P. velifera and P. gillii. Pipeline is implemented in snakemake form. 

To recreate, please obtain the reference genome at (INSERT RERFERENCE SOURCE HERE) and modify the data/refInfo.tsv file with the appropriate path info.

This pipeline is designed to be run in a computing environment with SLURM job submission, where each individual snakejob will be submitted as a separate SLURM job. The snakemake command should be run in a session that continues to operate in the background if the user is disconnected from the server (tmux is ideal for this). This is especially important as the variant calling process will take several days. Please view the runSnakemake.sh file for some example tmux and snakemake commands.