#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o testCactus

#set paths
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/cactus-bin-v2.9.9
CACTUS_PATH=/ocean/projects/bio230047p/mchin/software/cactus-bin-v2.9.9

#load python env
source $CACTUS_PATH/venv-cactus-v2.9.9/bin/activate

#make directory in scriptOuts for job outputs
mkdir scriptOuts/testMammalsSLURM

#run Cactus
cactus data/testRunCactus $CACTUS_PATH/examples/evolverMammals.txt outputs/testMammals/testMammalsSLURM.hal \
	--batchSystem slurm \
	--consCores 128 \
	--doubleMem true \
	--batchLogsDir scriptOuts/testMammalsSLURM \
	--workDir data/tmp \
	--maxMemory 256GB \
	--slurmTime 1:00:00 \
	--slurmPartition RM \
	--slurmGPUPartition RM-512

