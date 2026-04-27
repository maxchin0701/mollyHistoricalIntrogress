#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#pop0=Plat
#pop1=Pfor
#pop2=Pmex

#get model
model=(${snakemake_wildcards[modelAM]})

#get into right wd
cd data/${model}

#symlink
ln -s ../sfs/demModelMolly_jointMAFpop1_0.obs ./${model}_jointMAFpop1_0.obs
ln -s ../sfs/demModelMolly_jointMAFpop2_0.obs ./${model}_jointMAFpop2_0.obs
ln -s ../sfs/demModelMolly_jointMAFpop2_1.obs ./${model}_jointMAFpop2_1.obs
ln -s ../sfs/demModelMolly_jointMAFpop3_0.obs ./${model}_jointMAFpop3_0.obs
ln -s ../sfs/demModelMolly_jointMAFpop3_1.obs ./${model}_jointMAFpop3_1.obs
ln -s ../sfs/demModelMolly_jointMAFpop3_2.obs ./${model}_jointMAFpop3_2.obs






