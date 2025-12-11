#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o testHal2MAF

#set paths
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/cactus-bin-v2.9.9
CACTUS_PATH=/ocean/projects/bio230047p/mchin/software/cactus-bin-v2.9.9

#load python env
source $CACTUS_PATH/venv-cactus-v2.9.9/bin/activate

#run cactus_hal2maf
cactus-hal2maf data/testHAL2MAF outputs/testMammalsSLURM.hal outputs/testMammalsSLURM.maf \
	--refGenome simHuman_chr6 \
	--chunkSize 100000 \
	--batchCount 5 \
	--batchCores 128 \
	--filterGapCausingDupes \
	--dupeMode single \
	--batchSystem slurm \
	--batchLogsDir scriptOuts/testHAL2MAFSLURM \
	--slurmTime 1:00:00 \
	--slurmPartition RM \
	--maxMemory 256GB



