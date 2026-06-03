#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o scriptOuts/runSFSStatsSimExon

#initialize conda
module load anaconda3/2022.10
conda activate hts_popgen

#export to path
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/hts_popgen

for i in {1..1000}; do
	echo "Processing iteration ${i}"
	#variables
	POPS=data/popsHTSPopgen.tsv
	VCF=data/mergedHapCalls/parHapCalls.vcf.gz
	BED=output/LOHRegionsSharedSim/simLOHSharedExon_${i}Sorted.bed
	OUT=output/selectionStats/simStatsOut/simRegionsStatsExon_${i}.tsv

	#run SFS stats
	sfs --vcf ${VCF} --pops ${POPS} --bed ${BED} --outgroup Pret > ${OUT}
done