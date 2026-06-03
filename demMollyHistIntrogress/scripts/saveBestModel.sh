#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#pop0=Plat
#pop1=Pfor
#pop2=Pmex

#get model prefix
models=("noGeneFlow" "ancestralGeneFlow" "modernGeneFlow" "ancestralModernGeneFlowCons" "ancestralModernGeneFlowShift")

for model in ${models[@]}; do
	#get into model home dir
	cd data/${model}

	#retreieve iteration #, raw estimated lhood, and obs-est lhood
	bestrunLhoods=$(cat iter{1..100}/${model}/${model}.bestlhoods | grep -v MaxObsLhood | awk -v model="$model" '{print model,NR,$(NF-1),$NF-$(NF-1)}' | sort -k 4 | head -1)

	#add to file
	echo $bestrunLhoods >> ../modelBestLhoods.txt

	#get back into project root
	cd ../..
done

bestOverallLhood=$(cat data/modelBestLhoods.txt | sort -k 4 | head -1)

#save best model
bestModel=$(echo $bestOverallLhood | cut -f 1 -d " ")
#save best iteration
bestRun=$(echo $bestOverallLhood | cut -f 2 -d " ")
#save best lhood diff
bestLhoodDiff=$(echo $bestOverallLhood | cut -f 4 -d " ")

#copy best run to new dir
#mkdir bestIter
for model in ${models[@]}; do
	cp -r data/${model}/iter${bestRun}/* data/${model}/bestIter
	touch data/${model}/bestIter/${model}_temp.lhoods
done

#echo message
echo Best run: ${bestModel}, Iteration ${bestRun}, deltaLhood: ${bestLhoodDiff}
