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
cd ../data/rawSeq/"${snakemake_wildcards[sample]}"

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
	mv "${snakemake_wildcards[sample]}".sra "$curSrx"/"$curSrx".sra

	fasterq-dump "$curSrx" --include-technical -S --outdir ./"$curSrx" --outfile "${snakemake_wildcards[sample]}"
	echo "successful fasterq-dump"
	#remove unnecessary files
	echo "gzip ./$curSrx\/"${snakemake_wildcards[sample]}"*\.fastq now"
	gzip ./$curSrx\/"${snakemake_wildcards[sample]}"*\.fastq
	echo "move now"
	mv ./$curSrx\/"${snakemake_wildcards[sample]}"*\.gz .
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
			mv "${snakemake_wildcards[sample]}".sra "$curSrx"/"$curSrx".sra
		else
			mv "${snakemake_wildcards[sample]}"\_$i.sra "$curSrx"/"$curSrx".sra
		fi

		fasterq-dump "$curSrx" --include-technical -S --outdir ./"$curSrx" --outfile "${snakemake_wildcards[sample]}"

		mv $curSrx\/"${snakemake_wildcards[sample]}"_1.fastq ./"${snakemake_wildcards[sample]}"_$i\_1.fastq
		mv $curSrx\/"${snakemake_wildcards[sample]}"_2.fastq ./"${snakemake_wildcards[sample]}"_$i\_2.fastq

		rm -r $curSrx
	done

	#concat forward and reverse
	cat *_1.fastq > "${snakemake_wildcards[sample]}"\_1.fastq 
	rm -r "${snakemake_wildcards[sample]}"\_*\_1.fastq
	cat *_2.fastq > "${snakemake_wildcards[sample]}"\_2.fastq 
	rm -r "${snakemake_wildcards[sample]}"\_*\_2.fastq
	echo "gzip "${snakemake_wildcards[sample]}"*\.fastq now"
	#gzip
	gzip "${snakemake_wildcards[sample]}"*\.fastq
fi