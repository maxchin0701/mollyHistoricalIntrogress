#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#pop0=Plat
#pop1=Pfor
#pop2=Pmex

#export fsc path
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/fsc28_linux64

#load in variables
model=(${snakemake_wildcards[model]})

#get into right wd
cd data/${model}/bestIter

#copy sfs
cp ${model}_jointMAFpop1_0.obs ${model}_maxL_jointMAFpop1_0.obs
cp ${model}_jointMAFpop2_0.obs ${model}_maxL_jointMAFpop2_0.obs
cp ${model}_jointMAFpop2_1.obs ${model}_maxL_jointMAFpop2_1.obs
cp ${model}_jointMAFpop3_0.obs ${model}_maxL_jointMAFpop3_0.obs
cp ${model}_jointMAFpop3_1.obs ${model}_maxL_jointMAFpop3_1.obs
cp ${model}_jointMAFpop3_2.obs ${model}_maxL_jointMAFpop3_2.obs

#copy maxL
cp ${model}/${model}_maxL.par .

for i in {1..100}
do
	echo Running iteration $i
	#run fsc
	fsc28 -i ${model}_maxL.par -m -0 -q -n 1000000 --logprecision 18 -c4 -B 8

	#append to lhood file 
	sed -n '2,3p' ${model}_maxL/${model}_maxL.lhoods  >> ${model}_temp.lhoods

	rm -r ${model}_maxL
done

#move lhood file
mv ${model}_temp.lhoods ${model}.lhoods