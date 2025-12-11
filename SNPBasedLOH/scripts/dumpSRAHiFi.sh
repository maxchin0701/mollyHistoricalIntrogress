#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 4:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o dumpSRA

#activate conda env
module load anaconda3/2022.10
conda activate braker

#get into script wd
cd scripts

#export path to sratoolkit
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/sratoolkit.3.1.1-centos_linux64/bin

#extract reads
cd ../data/rawSeq/"${snakemake_wildcards[sample_hifi]}"

#assign snakemake srx param to var
srx=(${snakemake_params[srx]})

echo "${#srx[@]}"

#check if multiple sra files
if [ "${#srx[@]}" == "1" ] ;
then
	#set curSRX
	curSrx=${srx[0]}

	#mv to new directory
	mkdir "$curSrx"
	mv "${snakemake_wildcards[sample_hifi]}"_HiFi.sra "$curSrx"/"$curSrx".sra

	fasterq-dump "$curSrx" --outdir ./"$curSrx" --outfile "${snakemake_wildcards[sample_hifi]}".fastq --threads 10
	echo "successful fasterq-dump"
	#remove unnecessary files
	echo "gzip ./$curSrx\/"${snakemake_wildcards[sample_hifi]}".fastq now"
	gzip ./$curSrx\/"${snakemake_wildcards[sample_hifi]}".fastq
	echo "move now"
	mv ./$curSrx\/"${snakemake_wildcards[sample_hifi]}"*\.gz ./"${snakemake_wildcards[sample_hifi]}"_Hifi.fastq.gz
	rm -r "$curSrx"
else
	for (( i=0; i<${#srx[@]}; i++ ));
	do
		#set curSRX
		curSrx=${srx[$i]}

		#mv to new directory
		mkdir "$curSrx"
		if [ "$i" == "0" ] ;
		then
			mv "${snakemake_wildcards[sample_hifi]}"_.sra "$curSrx"/"$curSrx".sra
		else
			mv "${snakemake_wildcards[sample_hifi]}"\_$i.sra "$curSrx"/"$curSrx".sra
		fi

		fasterq-dump "$curSrx" --outdir ./"$curSrx" --outfile "${snakemake_wildcards[sample_hifi]}"

		mv $curSrx\/"${snakemake_wildcards[sample_hifi]}".fastq ./"${snakemake_wildcards[sample_hifi]}"_$i\_Hifi.fastq

		rm -r $curSrx
	done

	#concat forward and reverse
	cat *_Hifi.fastq > "${snakemake_wildcards[sample_hifi]}"\_Hifi.fastq 
	rm -r "${snakemake_wildcards[sample_hifi]}"\_*\_Hifi.fastq
	echo "gzip "${snakemake_wildcards[sample_hifi]}"_Hifi.fastq now"
	#gzip
	gzip "${snakemake_wildcards[sample_hifi]}"_Hifi.fastq
fi