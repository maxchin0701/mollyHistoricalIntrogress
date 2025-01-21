#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o createGenomicsDB

#activate conda env
module load GATK
module load anaconda3/2022.10
conda activate genomics

#collect sample vcf files
gatk GenomicsDBImport \
	--genomicsdb-workspace-path data/hapCalls/genDB_"${snakemake_wildcards[chr]}" \
	--sample-name-map "${snakemake_input[0]}" \
	-L "${snakemake_input[1]}"

#touch output file
touch ${snakemake_output[0]}
