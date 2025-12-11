#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#pop0=Plat
#pop1=Pfor
#pop2=Pmex

#activate conda env
module load GATK
module load anaconda3/2022.10
conda activate genomics

#set bam list
vcfList=(${snakemake_input[0]})
seqDict=(${snakemake_input[1]})
samps=(${snakemake_input[2]})

#run mergeVCF from gatk
gatk MergeVcfs -I $vcfList \
	-D $seqDict \
	-O data/combinedHap/chrAllCombinedFilt.vcf.gz

#subset samples from vcf
gatk SelectVariants -V data/combinedHap/chrAllCombinedFilt.vcf.gz \
	-sn $samps \
	-O data/combinedHap/chrSymParCombinedFilt.vcf.gz
