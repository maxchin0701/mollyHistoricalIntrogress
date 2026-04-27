#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o subsetVCFParHyb

#pop0=Plat
#pop1=Pfor
#pop2=Pmex

#initialize conda
source ~/.bashrc

#activate conda env
module load gatk/4.5.0.0
module load bcftools/1.19
conda activate gatk

#get samps
SAMPS=data/sampsParHyb.args
CHR="${snakemake_wildcards[chr]}"

#check for gzip
if ! file "data/hapCalls/${CHR}CombinedFiltParHyb.vcf.gz" | grep -q "gzip compressed"; then
    mv data/hapCalls/${CHR}CombinedFiltParHyb.vcf.gz data/hapCalls/${CHR}CombinedParHyb.vcf
    bgzip data/hapCalls/${CHR}CombinedFiltParHyb.vcf
fi

#index vcf file
tabix -p vcf data/hapCalls/${CHR}CombinedParHyb.vcf.gz

#subset samples from vcf
gatk SelectVariants -V data/hapCalls/${CHR}CombinedParHyb.vcf.gz \
	-sn $SAMPS \
	-O data/hapCallParHybs/${CHR}Unfilt.vcf.gz

#filter VCF
bcftools filter \
	data/hapCallParHybs/${CHR}Unfilt.vcf.gz \
	-i 'QUAL>30 && Type="snp"' | \
	bcftools view -c1 -m2 -v snps -Oz -o data/hapCallParHybs/${CHR}.vcf.gz

tabix -p vcf data/hapCallParHybs/${CHR}.vcf.gz