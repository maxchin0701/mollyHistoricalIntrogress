#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 4:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o downloadSRA

#activate conda env
module load anaconda3/2022.10
conda activate braker

#get into scripts directory
cd scripts

#export path to sratoolkit
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/sratoolkit.3.1.1-centos_linux64/bin

#assign snakemake srx param to var
srx=(${snakemake_params[srx]})

#check if multiple sra files

if [ "${#srx[@]}" == "1" ] ;
then
	#set curSRX
	curSrx=${srx[0]}

	#echo
	echo "Downloading $curSrx"

	#prefetch reads
	prefetch "$curSrx" -O ../data/rawSeq

	#mv to new directory
	cd ../data/rawSeq
	mv "$curSrx"/* "${snakemake_wildcards[sample]}"
	mv "${snakemake_wildcards[sample]}"/"$curSrx".sra "${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}".sra
	rm -r "$curSrx"
else 
	for (( i=0; i<${#srx[@]}; i++ ));
	do
		curSrx=${srx[$i]}

		#echo
		echo "Downloading $curSrx"

		#prefetch reads
		prefetch "$curSrx" -O ../data/rawSeq

		#mv to new directory
		cd ../data/rawSeq
		mv "$curSrx"/* "${snakemake_wildcards[sample]}"

		#rename appropriately
		if [ "$i" == "0" ] ;
		then
			mv "${snakemake_wildcards[sample]}"/"$curSrx".sra "${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}".sra
		else
			mv "${snakemake_wildcards[sample]}"/"$curSrx".sra "${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\_$i.sra
		fi

		rm -r "$curSrx"

		cd ../../scripts
	done
fi
