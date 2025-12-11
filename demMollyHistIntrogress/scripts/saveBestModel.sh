#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#pop0=Plat
#pop1=Pfor
#pop2=Pmex

#get model prefic
model=(${snakemake_wildcards[model]})

#get into model home dir
cd data/${model}

#retreieve iteration #, raw estimated lhood, and obs-est lhood
bestrunLhoods=$(cat iter{1..100}/${model}/${model}.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$(NF-1),$NF-$(NF-1)}' | sort -k 2 | head -1)

#save best iteration
bestRun=$(echo $bestrunLhoods | cut -f 1 -d " ")
#save best lhood
bestLhood=$(echo $bestrunLhoods | cut -f 2 -d " ")
#save best lhood diff
bestLhoodDiff=$(echo $bestrunLhoods | cut -f 3 -d " ")

#copy best run to new dir
#mkdir bestIter
cp -r iter${bestRun}/* bestIter
touch bestIter/${model}_temp.lhoods

#echo message
echo Best run: Iteration ${bestRun}
