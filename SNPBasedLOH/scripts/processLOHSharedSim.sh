#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 30:00:00
#SBATCH --ntasks-per-node=36
#SBATCH -o align

#load modules
module load bedtools

#get into right directory
cd output/LOHRegionsSharedSim

for i in  $(ls -f *.*); do
(

	prefix=$(echo $i | cut -d '.' -f 1 )
	#sort regions
	sort -k1,1 -k2,2n $prefix\.bed > $prefix\Sorted.bed

	#bedtools intersect
	bedtools intersect -a $prefix\Sorted.bed \
		-b ../../data/funcAnnotate/eMapperOut/funcAnnotatePFor.emapper.decorated.gff -wa -wb -sorted \
		-g ../../data/PForHapMex.genome > ../LOHRegionsSimIntersectOut/$prefix\Out.tsv
)
done