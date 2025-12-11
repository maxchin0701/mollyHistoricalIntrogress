#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 10:00:00
#SBATCH --ntasks-per-node=30
#SBATCH -o ploidyAll

#activate conda env
module load anaconda3/2022.10
conda activate gatk

#path for local gatk
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/gatk-4.1.9.0/

#take in snakemake vars
ref="${snakemake_input[ref]}"
readCounts="${snakemake_input[readCounts]}"
intervals="${snakemake_input[intervals]}"
outFile="${snakemake_output[0]}"
chr="${snakemake_wildcards[chr]}"

#create ploidy prior files
python scripts/createPloidyProb.py --inFasta $ref --scaff $chr --dipProb 1 --tripProb 0 --outTSV data/gCNV/ploidy/scaffPloidyPriors_$chr\.tsv

#build list of input samp read counts
readCountsAll=""

for sampReadCounts in $readCounts; do 
	readCountsAll+="--input $sampReadCounts "
done

#determine contig ploidy
gatk DetermineGermlineContigPloidy \
	$readCountsAll\
	--contig-ploidy-priors data/gCNV/ploidy/scaffPloidyPriors_$chr\.tsv \
	-L $intervals \
  	--interval-merging-rule OVERLAPPING_ONLY \
  	--output-prefix ploidy_$chr \
  	--output data/gCNV/ploidy


#create file if directoral exists
if [ -d "data/gCNV/ploidy/ploidy_$chr-calls" ]; then
	touch $outFile
fi

