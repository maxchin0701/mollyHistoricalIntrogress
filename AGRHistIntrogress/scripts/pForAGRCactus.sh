#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=
#SBATCH -o testCactus

#set paths
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/cactus-bin-v2.9.9
CACTUS_PATH=/ocean/projects/bio230047p/mchin/software/cactus-bin-v2.9.9

#load python env
source $CACTUS_PATH/venv-cactus-v2.9.9/bin/activate

#make directory in scriptOuts for job outputs
mkdir scriptOuts/pForAGRCactus

#run Cactus
cactus data/pForAGRCactus data/pForAGRCactusConfig/pForAGRCactusConfig.txt outputs/pForAGRCactus/pForAGRCactus.hal \
	--batchSystem slurm \
	--consCores 64 \
	--doubleMem true \
	--batchLogsDir scriptOuts/pForAGRCactus \
	--workDir data/tmp \
	--maxMemory 128GB \
	--slurmTime 24:00:00 \
	--slurmPartition RM \
	--slurmGPUPartition RM-512

