#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o scriptOuts/generateTinkerRefPanel

#initialize conda
source ~/.bashrc

#load environment
module load bcftools/1.19
module load tabix/2013-12-16
conda activate R

#set chromosome
CHR="${snakemake_wildcards[chr]}"

#pseudophase with RScript
Rscript scripts/pseudoPhaseHyb.R ${CHR}

#unzip reference panels
gunzip data/refPanel/${CHR}WhatshapRefPanel.vcf.gz

#bgzip reference panels
bgzip data/refPanel/${CHR}WhatshapRefPanel.vcf

#index vcf reference panels
bcftools tabix data/refPanel/${CHR}WhatshapRefPanel.vcf.gz


