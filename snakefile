#import pandas df
import pandas as pd

#read in samples, reference, chromosomes
#samp_df = pd.read_csv("./data/testSampleIndex.tsv", sep='\t') 
samp_df = pd.read_csv("./data/sampleIndex.tsv", sep='\t') 
ref_df = pd.read_csv("./data/refInfo.tsv", sep="\t")
chr_df = pd.read_csv("./data/chrIndex.tsv", sep="\t")
filtCut_df = pd.read_csv("./data/filterCutoffs.tsv", sep="\t")
refSpecies = ["LM","GLMV"]
#chr_df = pd.read_csv("./data/testChrIndex.tsv", sep="\t")
#filtCut_df = pd.read_csv("./data/testFilterCutoffs.tsv", sep="\t")

#creating new df for combined cutoffs + chrs
#initialize vectors to store
filtCutPar=[]
filtCutChrs=[]
filtCutGQ=[]
filtCutDP=[]

#create separate group of rows for each 
for i in refSpecies:
    for j in chr_df["CHR"]:
        for k in range(0,len(filtCut_df)):
            filtCutPar.append(i)
            filtCutChrs.append(j)
            filtCutGQ.append(filtCut_df["GQ"][k])
            filtCutDP.append(filtCut_df["DP"][k])
    
#append to create new df
filtCutFinal_df=pd.DataFrame({"PARENTS": filtCutPar, "CHR": filtCutChrs, "GQ": filtCutGQ, "DP": filtCutDP})
    

rule all:
    input:
        expand("output/LOHRegions/{chr}LOHRegions_GQ{gqcut}_DP{dpcut}_{pars}.tsv",zip,pars = filtCutFinal_df["PARENTS"],chr = filtCutFinal_df["CHR"],gqcut = filtCutFinal_df["GQ"],dpcut = filtCutFinal_df["DP"]),
        expand("output/LOHAnc/{chr}LOHAnc_GQ{gqcut}_DP{dpcut}_{pars}.tsv",zip, pars = filtCutFinal_df["PARENTS"], chr = filtCutFinal_df["CHR"],gqcut = filtCutFinal_df["GQ"],dpcut = filtCutFinal_df["DP"]),
        expand("output/SNPall/{chr}SNPall_GQ{gqcut}_DP{dpcut}_{pars}.tsv",zip, pars = filtCutFinal_df["PARENTS"], chr = filtCutFinal_df["CHR"],gqcut = filtCutFinal_df["GQ"],dpcut = filtCutFinal_df["DP"]),
        expand("output/fixDiffVCF/{chr}FixDiff_GQ{gqcut}_DP{dpcut}_{pars}.vcf.gz", zip, pars = filtCutFinal_df["PARENTS"], chr = filtCutFinal_df["CHR"],gqcut = filtCutFinal_df["GQ"],dpcut = filtCutFinal_df["DP"]),
        expand("output/haplostrips/haplostrips_{chr}_GQ{gqcut}_DP{dpcut}_{pars}.pdf",zip, pars = filtCutFinal_df["PARENTS"], chr = filtCutFinal_df["CHR"],gqcut = filtCutFinal_df["GQ"],dpcut = filtCutFinal_df["DP"]),


#download sra files
rule prefetchSRA:
    priority: 10
    params:
        srx=lambda wc: (samp_df[samp_df['NAME'] == wc.sample]['SRX']).tolist(),
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=60,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=10 -o scriptOuts/prefetchSRA_{wildcards.sample} -t 1:00:00",
    output:
        "data/rawSeq/{sample}/{sample}.sra",
    threads: 10
    log: 
        "data/rawSeq/{sample}/prefetchSRA_{sample}.log"
    script:
        "scripts/prefetchSRA.sh"

#extract fastq.gz from sra files
rule dumpSRA:
    priority: 9
    input:
        "data/rawSeq/{sample}/{sample}.sra",
    params:
        srx=lambda wc: (samp_df[samp_df['NAME'] == wc.sample]['SRX']).tolist(),
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=240,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=10 -o scriptOuts/dumpSRA_{wildcards.sample} -t 4:00:00",
    output:
        "data/rawSeq/{sample}/{sample}_1.fastq.gz",
        "data/rawSeq/{sample}/{sample}_2.fastq.gz",
    threads: 10
    log: 
        "data/rawSeq/{sample}/dumpSRA_{sample}.log"
    script:
        "scripts/dumpSRA.sh"
        
#trim files
rule trimFastq:
    priority: 8
    input:
        forRead="data/rawSeq/{sample}/{sample}_1.fastq.gz",
        revRead="data/rawSeq/{sample}/{sample}_2.fastq.gz",
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=90,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=4 -o scriptOuts/trimFastq_{wildcards.sample} -t 1:30:00",
    output:
        "data/trimSeq/{sample}/{sample}_forward_paired.fastq.gz",
        "data/trimSeq/{sample}/{sample}_reverse_paired.fastq.gz",
    threads: 4
    log: 
        "data/trimSeq/{sample}/trimFastq_{sample}.log"
    script:
        "scripts/trimFastq.sh"
 
#index ref
rule indexRef:
    priority: 10
    input:
        "data/refGenome/{refName}.fasta",
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=60,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=20 -o scriptOuts/indexRef -t 1:00:00",    
    output:
        "data/refGenome/{refName}.fasta.fai",
        "data/refGenome/{refName}.dict"
    threads: 20
    log: 
        "data/refGenome/indexRef_{refName}.log"
    script:
        "scripts/indexRef.sh"

#create interval lists for chrs
rule createIntervalL:
    priority: 9
    input:
        "data/refGenome/PForSch01.fasta",
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=60,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=1 -o scriptOuts/createIntervalL_{wildcards.chr} -t 1:00:00",    
    output:
        "data/refGenome/intervalLists/PFor_{chr}.bed",
        "data/refGenome/intervalLists/PFor_{chr}.interval_list"
    threads: 1
    log: 
        "data/refGenome/intervalLists/createIntervalL_{chr}.log"
    script:
        "scripts/createIntervalL.sh"
        
#create genome wide interval list
rule createGWideIntervalL:
    priority: 8
    input:
        bedF=expand("data/refGenome/intervalLists/PFor_{chr}.bed",chr = chr_df["CHR"])
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=60,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=1 -o scriptOuts/q -t 1:00:00",    
    output:
        "data/refGenome/intervalLists/PFor_genome.bed",
        "data/refGenome/intervalLists/PFor_genome.interval_list",
    threads: 1
    log: 
        "data/refGenome/intervalLists/createGWideIntervalL.log"
    script:
        "scripts/createGWideIntervalL.sh"

#run alignments
rule alignFastq:
    priority: 7
    input:
        forRead="data/trimSeq/{sample}/{sample}_forward_paired.fastq.gz",
        revRead="data/trimSeq/{sample}/{sample}_reverse_paired.fastq.gz",
        ref=ref_df.iloc[0,1],
        refIndex=expand("data/refGenome/{refName}.fasta.fai",refName = ref_df["NAME"])
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=240,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=36 -o scriptOuts/align_{wildcards.sample} -t 4:00:00",
    output:
        "data/aligned/{sample}/{sample}SortedDupMarked.bam"
    threads: 36
    log: 
        "data/aligned/{sample}/alignFastq_{sample}.log"
    script:
        "scripts/alignFastq.sh"

#call variants
rule callHaps:
    priority: 6
    input:
        bam="data/aligned/{sample}/{sample}SortedDupMarked.bam",
        intervals="data/refGenome/intervalLists/PFor_{chr}.interval_list",
        ref=ref_df.iloc[0,1],
        refIndex=expand("data/refGenome/{refName}.fasta.fai",refName = ref_df["NAME"])
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=600,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=10 -o scriptOuts/callHaps_{wildcards.sample}_{wildcards.chr} -t 10:00:00",
    output:
        "data/hapCalls/{sample}/{sample}{chr}.vcf.gz"
    threads: 10
    log: 
        "data/hapCalls/{sample}/callHaps_{sample}_{chr}.log"
    script:
        "scripts/callHaplotypes.sh"

#write sample_map 
rule createSampleMap:
    priority: 5
    input:
        sampVCF=expand("data/hapCalls/{sample}/{sample}{{chr}}.vcf.gz",sample = samp_df["NAME"].unique())
    params:
        sampNames=samp_df["NAME"].unique()
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=60,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=1 -o scriptOuts/createSampleMap_{wildcards.chr} -t 1:00:00",
    output:
        "data/hapCalls/allSamps_{chr}.sample_map"
    threads: 1
    log: 
        "data/hapCalls/createSampleMap_{chr}.log"
    script:
        "scripts/createSampleMap.sh"

#create genomicsDB database 
rule createGenomicsDB:
    priority: 4
    input:
        "data/hapCalls/allSamps_{chr}.sample_map",
        "data/refGenome/intervalLists/PFor_{chr}.interval_list",
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=240,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=10 -o scriptOuts/createGenomicsDB_{wildcards.chr} -t 4:00:00",
    output:
        "data/hapCalls/genDBAllExists_{chr}.txt"
    threads: 10
    log: 
        "data/hapCalls/createGenomicsDB_{chr}.log"
    script:
        "scripts/createGenomicsDB.sh"

#perform joint genotyping
rule jointGenotyping:
    priority: 3
    input:
        "data/hapCalls/genDBAllExists_{chr}.txt",
        ref=ref_df.iloc[0,1],
    params:
        genDBDir=lambda wildcards: "data/hapCalls/genDB_" + wildcards.chr,
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=120,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=10 -o scriptOuts/jointGenotyping_{wildcards.chr} -t 2:00:00",
    output:
        "data/hapCalls/{chr}Combined.vcf.gz",
        "data/hapCalls/{chr}CombinedFilt.vcf.gz"
    threads: 10
    log:
        "data/hapCalls/jointGenotyping_{chr}.log"
    script:
        "scripts/jointGenotyping.sh"

#run Rscript to collect LOH regions
rule collectLOH:
    priority: 2
    input:
        "data/hapCalls/{chr}CombinedFilt.vcf.gz"
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=60,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=10 -o scriptOuts/collectLOH_{wildcards.chr}_GQ{wildcards.gqcut}_DP{wildcards.dpcut}_{wildcards.pars} -t 1:00:00",
    output:
        "output/LOHRegions/{chr}LOHRegions_GQ{gqcut}_DP{dpcut}_{pars}.tsv",
        "output/LOHAnc/{chr}LOHAnc_GQ{gqcut}_DP{dpcut}_{pars}.tsv",
        "output/SNPall/{chr}SNPall_GQ{gqcut}_DP{dpcut}_{pars}.tsv",
        "output/fixDiffVCF/{chr}FixDiff_GQ{gqcut}_DP{dpcut}_{pars}.vcf.gz"
    threads: 10
    log:
        "output/collectLOH_{chr}_GQ{gqcut}_DP{dpcut}_{pars}.log"
    script:
        "scripts/findLOH.R"

#generate haplostrips
rule genHaplostrips:
    priority: 1
    input:
        "output/fixDiffVCF/{chr}FixDiff_GQ{gqcut}_DP{dpcut}_{pars}.vcf.gz"
    params:
        chrL=lambda wildcards: str(chr_df[chr_df['CHR'] == wildcards.chr]['END'].iloc[0]),
    resources:
        slurm_partition="RM-shared",
        nodes=1,
        runtime=60,
        slurm_extra=lambda wildcards: f"--ntasks-per-node=10 -o scriptOuts/genHaplostrips_{wildcards.chr}_GQ{wildcards.gqcut}_DP{wildcards.dpcut}_{wildcards.pars} -t 1:00:00",
    output:
        "output/haplostrips/haplostrips_{chr}_GQ{gqcut}_DP{dpcut}_{pars}.pdf",
    threads: 10
    log:
        "output/haplostrips/genHaplostrips_{chr}_GQ{gqcut}_DP{dpcut}_{pars}.log"
    script:
        "scripts/genHaplostrips.sh"



    




