#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 10:00:00
#SBATCH --ntasks-per-node=30
#SBATCH -o gCNVAll

#activate conda env
module load anaconda3/2022.10
conda activate gatk

#path for local gatk
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/gatk-4.1.9.0/

#take in snakemake vars
readCounts="${snakemake_input[readCounts]}"
intervals="${snakemake_input[intervals]}"
outFile="${snakemake_output[0]}"
chr="${snakemake_wildcards[chr]}"

#build list of input samp read counts
readCountsAll=""

for sampReadCounts in $readCounts; do 
	readCountsAll+="--input $sampReadCounts "
done

#run intervals
gatk GermlineCNVCaller \
	$readCountsAll \
	--run-mode COHORT \
	-L $intervals \
	--interval-merging-rule OVERLAPPING_ONLY \
	--contig-ploidy-calls data/gCNV/ploidy/ploidy_$chr\-calls \
	--output data/gCNV/gCNV \
	--output-prefix gCNV_$chr

#create file if directoral exists
if [ -d "data/gCNV/gCNV/gCNV_$chr-calls" ]; then
	touch $outFile
fi
