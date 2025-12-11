#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 11:33:18 2025

@author: maxchin
"""

#import modules
from Bio import AlignIO
from Bio import SeqIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
import argparse
import pandas as pd

#function for getting postion in alignment
def getAlleleAtPos(alignment,ungapLength,refSampIndex):
    #get reference sequence
    refSeq=alignment[refSampIndex].seq
    #get gapped length
    curUngapLength=0
    curGapLength=0
    while curUngapLength < ungapLength:
        if refSeq[curGapLength] != "-":
            #print(refSeq[curGapLength])
            curUngapLength+=1
            #print(curUngapLength)
        curGapLength+=1
        #print(curGapLength,"\n")
    
    #assign final ungap index
    snpGapIndex=curGapLength
    
    #get bases at position for each sample
    alleles=[]
    for sampSeq in curAlign:
        alleles.append(sampSeq.seq[snpGapIndex].capitalize())
        
    #return alleles
    return(alleles)

#parse user inputs
parser = argparse.ArgumentParser(description="Parsing test MAF output")
parser.add_argument("--inMAF", action="store", dest="inMAF", type=str)
#parser.add_argument("--inFASTA", action="store", dest="inFASTA", type=str)
parser.add_argument("--inIndex", action="store", dest="inIndex", type=str)
parser.add_argument("--indexSamp", action="store", dest="indexSamp", type=str)
parser.add_argument("--inVCF", action="store", dest="inVCF", type=str)
parser.add_argument("--inBED", action="store", dest="inBED", type=str)
parser.add_argument("--inChrs", action="store", dest="inChrs", type=str)
args=parser.parse_args()


#read in chromosomes
chrDf=pd.read_csv(args.inChrs, sep="\t")
chroms=chrDf.iloc[:,0]
#chroms=["chr1"]

#objects to store length
lenUngap=[]
lenGap=[]
"""
#loop through chromosomes
for chrom in chroms:
    #print chromosome
    print("traversing",chrom)
        
    #load in current chromosomes 
    curChromAlignments=list(AlignIO.parse(args.inMAF+"/"+chrom+".maf", "maf"))
    
    #get average length of alignments
    for curAlign in curChromAlignments:        
        lenUngap.append(curAlign[0].annotations["size"])
        lenGap.append(len(curAlign[0].seq))
        
print("Average alignment ungapped length",sum(lenUngap)/len(lenUngap),"total alignments:",len(lenUngap))
print("Average alignment gapped length",sum(lenGap)/len(lenGap))
"""

#read in bed and  files
LOHTracts = pd.read_csv(args.inBED, sep="\t",header=None)
LOHTest = LOHTracts[LOHTracts.iloc[:,0]=="chr1"]
LOHSnps = pd.read_csv(args.inVCF, sep="\t",header=None)

#counter variables
nCompleteRatio=0
nAGRSupport=0

#loop through LOH shared tracts
for i in range(len(LOHTracts)):
    #get current chromosome, tract start and end
    chrom=LOHTracts.iloc[i,0]
    tractStart=LOHTracts.iloc[i,1]
    tractEnd=LOHTracts.iloc[i,2]
    
    #subset SNPs
    curTractSnps=LOHSnps[(LOHSnps.iloc[:,4] == chrom) & 
                         (LOHSnps.iloc[:,5] >= tractStart) & 
                         (LOHSnps.iloc[:,5] <= tractEnd)]
    
    #read in MAF index
    curIndex = MafIO.MafIndex((args.inIndex + "/" + chrom + ".mafindex"),
                               (args.inMAF + "/" + chrom + ".maf"),
                               (args.indexSamp + "." + chrom))
    
    fixDiffSnps=len(curTractSnps)
    fixDiffAGR=0
    LOHAncestral=0
    #loop through SNPs
    for j in range(len(curTractSnps)):
        #print("\nnew SNP")
        curOverlaps=curIndex.search([curTractSnps.iloc[j,5]-1],
                                     [curTractSnps.iloc[j,5]])
        for curAlign in curOverlaps:
            if(len(curAlign) != 9): 
                continue
            #retreive samp idxs
            hapMexIdx=[k for k, curSamp in enumerate(curAlign) if 
                          curSamp.id.split(".")[0] == 'PForHapMex'][0]
            hapLatIdx=[k for k, curSamp in enumerate(curAlign) if 
                          curSamp.id.split(".")[0] == 'PForHaplat'][0]
            pMexIdx=[k for k, curSamp in enumerate(curAlign) if 
                          curSamp.id.split(".")[0] == 'PMex'][0]
            pLatIdx=[k for k, curSamp in enumerate(curAlign) if 
                          curSamp.id.split(".")[0] == 'PLat'][0]
            ancestralHapMexIdx=[k for k, curSamp in enumerate(curAlign) if 
                          curSamp.id.split(".")[0] == 'AncPForHapMex'][0]
            ancestralHapLatIdx=[k for k, curSamp in enumerate(curAlign) if 
                          curSamp.id.split(".")[0] == 'AncPForHapLat'][0]

            #get snp index
            snpIdx=curTractSnps.iloc[j,5]-1
            
            #diff in index
            #print(curAlign)
            #print(curAlign[hapMexIdx].annotations["start"])
            print(snpIdx)
            #print(snpIdx-curAlign[hapMexIdx].annotations["start"])
            locAlleles=getAlleleAtPos(curAlign,
                         snpIdx-curAlign[hapMexIdx].annotations["start"],
                         hapMexIdx)
            print(locAlleles)

            if(locAlleles[pMexIdx] != locAlleles[pLatIdx] and
               locAlleles[hapMexIdx] == locAlleles[hapLatIdx]):
                #print("locus has fix diff")
                fixDiffAGR+=1

                if locAlleles[ancestralHapMexIdx] == locAlleles[ancestralHapLatIdx]:
                    LOHAncestral+=1
    
    if fixDiffAGR==0:
       print(chrom," ",tractStart,"-",tractEnd,": ",
             "No fixed differences found in tract",sep="")
    else :
        #add to counters
        if(fixDiffAGR/fixDiffSnps == 1):
            nCompleteRatio+=1
        
        if(LOHAncestral/fixDiffAGR == 1):
            nAGRSupport+=1
        
        #print statements
        print("Number of SNPs in BED: ",fixDiffSnps,"       Number of SNPs in AGR: ", fixDiffAGR, "       Ratio: ",fixDiffAGR/fixDiffSnps)
        print("Proportion of SNPs LOH in ancestral for tract ",chrom," ",tractStart,"-",tractEnd,": ",
              LOHAncestral/fixDiffAGR,sep="")

print("Number of tracts with all SNPs represented in AGR: ",nCompleteRatio/len(LOHTracts))
print("Number of tracts supported by AGR: ",nAGRSupport/len(LOHTracts))

                             
#Identify unique 
"""
#set chr
chroms=["chr1","chr2","chr3","chr4","chr5","chr6","chr7_24","chr8","chr9","chr10",
        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chr20","chr21","chr22","chr23"]
for i in chroms:
    #get current chroms
    chrom=i
    
    #build index
    testIndex = MafIO.MafIndex((args.inIndex + "/" + chrom + ".mafindex"),
                               (args.inMAF + "/" + chrom + ".maf"),
                               (args.indexSamp + "." + chrom))
    
    #counter variables
    numIncomplete=0
    numAlignments=0
    #sum up incomplete total alignments
    for align in AlignIO.parse((args.inMAF + "/" + chrom + ".maf"), "maf"):
        numAlignments+=1
        curLengths=[]
        if len(align) != 9:
            numIncomplete+=1
    #print
    #print(numAlignments)
    print(chrom,": ",(numIncomplete/numAlignments)*100,"% of alignments are incomplete")
"""