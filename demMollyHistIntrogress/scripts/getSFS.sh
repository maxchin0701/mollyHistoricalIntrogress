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
source ~/.bashrc
conda activate python3

#get population
#pop=(${snakemake_wildcards[popCombo]})
inVCF=(${snakemake_input[inVCF]})
pops=(${snakemake_input[pops]})
iter=(${snakemake_wildcards[iter]})

#create sfs dir if it doesn't exist
if [ ! -d data/sfs/iter${iter} ]; then
	mkdir data/sfs/iter${iter}
fi

#get real sfs
python3.14 scripts/buildMixedSFS.py --inVCF ${inVCF} --inSampleMap ${pops} --pop1 pop1 --pop2 pop0 --outSFS ${snakemake_output[0]}
python3.14 scripts/buildMixedSFS.py --inVCF ${inVCF} --inSampleMap ${pops} --pop1 pop2 --pop2 pop0 --outSFS ${snakemake_output[1]}
python3.14 scripts/buildMixedSFS.py --inVCF ${inVCF} --inSampleMap ${pops} --pop1 pop2 --pop2 pop1 --outSFS ${snakemake_output[2]}
python3.14 scripts/buildMixedSFS.py --inVCF ${inVCF} --inSampleMap ${pops} --pop1 pop3 --pop2 pop0 --outSFS ${snakemake_output[3]}
python3.14 scripts/buildMixedSFS.py --inVCF ${inVCF} --inSampleMap ${pops} --pop1 pop3 --pop2 pop1 --outSFS ${snakemake_output[4]}
python3.14 scripts/buildMixedSFS.py --inVCF ${inVCF} --inSampleMap ${pops} --pop1 pop3 --pop2 pop2 --outSFS ${snakemake_output[5]}


#rename
#mv data/sfs/${popCombo}/jointMAF${popCombo}.sfs data/sfs/${popCombo}/jointMAF${popCombo}.obs



