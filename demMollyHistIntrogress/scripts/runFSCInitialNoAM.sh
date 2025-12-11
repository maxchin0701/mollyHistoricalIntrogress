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
iter=(${snakemake_wildcards[iter]})
model=(${snakemake_wildcards[modelNoAM]})
#iter=1
#model=noGeneFlow

#get into right directory and make directory for iteration, clone into iter directory 
cd data/${model}
#mkdir iter${iter}
cd iter${iter}
if [ ! -f ${model}.est ]; then
	ln -s ../${model}* .
fi

#run fsc
fsc28 -t ${model}.tpl -e ${model}.est -m -0 -C 10 -n 1000000 -L 100 -s 0 -M --logprecision 18 -c4 -B 8 -y 3 -q

