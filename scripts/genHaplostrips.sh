#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o jointGenotyping

#activate conda env
module load anaconda3/2022.10
conda activate haplostrips

#export haplostrips path
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/haplostrips/

#change into dir
cd output/haplostrips

chr="${snakemake_wildcards[chr]}"
chrL="${snakemake_params[chrL]}"
gq="${snakemake_wildcards[gqcut]}"
dp="${snakemake_wildcards[dpcut]}"

#run haplostrips
haplostrips -v ../fixDiffVCF/$chr\FixDiff_GQ$gq\_DP$dp\.vcf.gz -i $chr\:1-$chrL -P ../../data/haplostripsPops.tsv -o haplostrips_$chr\_GQ$gq\_DP$dp -p plat,pmex,pfor -S 3 -C '#51127c','#fc8961','#b73779'