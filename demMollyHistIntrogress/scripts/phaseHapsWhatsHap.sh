#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 4:00:00
#SBATCH --ntasks-per-node=8
#SBATCH -o scriptOuts/phaseHapsWhatsHap

#initialize conda
source ~/.bashrc

#load conda env
module load gatk/4.5.0.0
module load bcftools/1.19
conda activate whatshap

#set chr
CHR="${snakemake_wildcards[chr]}"
SAMP="${snakemake_wildcards[samp]}"
REF="${snakemake_params[ref]}"
BAM=data/alignedPopSamples/${SAMP}/${SAMP}SortedDupMarked.bam
SCAFF="${snakemake_input[1]}"
VCF="${snakemake_input[0]}"
SAMP_VCF=data/hapCallParHybs/${CHR}_${SAMP}.vcf.gz

#subset sample
gatk SelectVariants -V ${VCF} \
	-sn ${SAMP} \
	-O data/hapCallParHybs/${CHR}_${SAMP}.vcf.gz

#index vcf
tabix -p vcf data/hapCallParHybs/${CHR}_${SAMP}.vcf.gz

#phase chromosome vcf
whatshap phase -o data/phasedHapsWhatshap/${CHR}_${SAMP}Phased.vcf.gz \
	--reference=${REF} ${SAMP_VCF} ${BAM} ${SCAFF}

#index chromosome vcf 
bcftools tabix data/phasedHapsWhatshap/${CHR}_${SAMP}Phased.vcf.gz

#remove temporary samp vcf files
rm -r data/hapCallParHybs/${CHR}_${SAMP}.vcf.gz*








