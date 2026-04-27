#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 18:21:24 2026

@author: maxchin
"""

#import modules
import argparse
import pandas as pd
from cyvcf2 import VCF
import numpy as np
from math import comb

def getProjectionVals(pop0Indices,pop1Indices,vcfPath,minRetention):
    #get number of chromosomes per locus for each pop
    
    #read in vcf
    vcf=VCF(vcfPath)
    
    #objects to store
    nPop0All = []
    nPop1All = []
    
    #loop through variants
    for variant in vcf:

        #pbjects to increment
        nPop0Cur = 0
        nPop1Cur = 0
        
        #loop through samples
        for samp in pop0Indices:
            #get sample genotype
            curGeno = variant.genotypes[samp]
            
            #increment current chr count
            nPop0Cur += sum([1 for all in curGeno[0:(len(curGeno) - 1)] if all != -1])
        
        #loop through samples
        for samp in pop1Indices:
            #get sample genotype
            curGeno = variant.genotypes[samp]
            #print(curGeno)
            
            #increment current chr count
            nPop1Cur += sum([1 for all in curGeno[0:(len(curGeno) - 1)] if all != -1])
        
        #append chr counts to list
        nPop0All.append(nPop0Cur)
        nPop1All.append(nPop1Cur)
    
    #get range of possible projection values
    pop0Project = range(int(max(nPop0All)/2), max(nPop0All)+1)
    pop1Project = range(int(max(nPop1All)/2), max(nPop1All)+1)
    
    #object to store projection candidate loci returned
    projResults = []
    
    #loop through projection values
    for proj0 in pop0Project:
        for proj1 in pop1Project:
            #retreive idxs for pop0 > proj0
            pop0LocIdx = [locIdx for locIdx, loc in enumerate(nPop0All) if loc >= proj0]
            
            #retreive idxs for pop1 > proj1
            pop1LocIdx = [locIdx for locIdx, loc in enumerate(nPop1All) if loc >= proj1]
            
            #intersect the two lists
            propRetainedPop0 = len(list(set(pop0LocIdx).intersection(pop1LocIdx)))/len(pop0LocIdx)
            propRetainedPop1 = len(list(set(pop0LocIdx).intersection(pop1LocIdx)))/len(pop1LocIdx)
            
            #append to results
            projResults.append([proj0, proj1, (propRetainedPop0 + propRetainedPop1)/2])
    
    #get best projection values
    bestProj = max([proj for proj in projResults if proj[2] >= minRetention], key=lambda x: x[0] + x[1])
    
    #return best projections
    return(bestProj[0:2])
                    
def getAlleleCounts(pop0Indices,pop1Indices,gts,allele):
    
    #list to store
    minorAlleleCounts=[]
    
    #subset genotypes
    subsetGtsPop0 = [gts[Idx] for Idx in pop0Indices]
    
    #subset genotypes
    subsetGtsPop0 = [gt[0:(len(gt) - 1)] for gt in subsetGtsPop0]
    
    #unlist gts
    allelesPop0 = [allele for gt in subsetGtsPop0 for allele in gt]
    
    #count minor alleles
    minorAlleleCounts.append(allelesPop0.count(allele))
    
    #subset genotypes
    subsetGtsPop1 = [gts[Idx] for Idx in pop1Indices]
    
    #subset genotypes
    subsetGtsPop1 = [gt[0:(len(gt) - 1)] for gt in subsetGtsPop1]
    
    #unlist gts
    allelesPop1 = [allele for gt in subsetGtsPop1 for allele in gt]
    
    #count minor alleles
    minorAlleleCounts.append(allelesPop1.count(allele))
    
    return(minorAlleleCounts)

def projectCounts(empMAC,sampSize,projSampSize):
    
    #array for storing projection probs
    projProbs = np.zeros(projSampSize + 1)
    
    #loop through possible projection values
    for projMAC in range(0,projSampSize + 1):
        projProbs[projMAC] = (
                comb(empMAC, projMAC) *
                comb(sampSize - empMAC, projSampSize - projMAC) /
                comb(sampSize, projSampSize)
            )
        
    #return projection values
    return(projProbs)
           
#take user input
parser = argparse.ArgumentParser(description="Generating 2D folded SFS from mixed-ploidy VCF")
parser.add_argument("--inVCF",required=True, action="store", dest="inVCF", type=str)
parser.add_argument("--inSampleMap",required=True, action="store", dest="inSampleMap", type=str)
parser.add_argument("--pop1",required=True, action="store", dest="pop1", type=str)
parser.add_argument("--pop2",required=True, action="store", dest="pop2", type=str)
parser.add_argument("--outSFS",required=True, action="store", dest="outSFS", type=str)
args=parser.parse_args()

#read in sample map
sampMap = pd.read_csv(args.inSampleMap, sep="\t",dtype=str)
sampMap['POP'] = sampMap['POP'].str.strip()

#get samples for population
pop0Samps = sampMap.loc[sampMap['POP'] == args.pop1]['SAMP'].tolist()
pop1Samps = sampMap.loc[sampMap['POP'] == args.pop2]['SAMP'].tolist()

#load in vcf
vcf=VCF(args.inVCF)

#get sample indices
pop0Indices = [sampIdx for sampIdx, samp in enumerate(vcf.samples) if samp in pop0Samps]
pop1Indices = [sampIdx for sampIdx, samp in enumerate(vcf.samples) if samp in pop1Samps]

#get projection values
projVals = getProjectionVals(pop0Indices,pop1Indices,args.inVCF,0.8)
pop0Proj = projVals[0]
pop1Proj = projVals[1]

#print projection values
print("projection value for population " + args.pop1 + " is " + str(pop0Proj))
print("projection value for population " + args.pop2 + " is " + str(pop1Proj))

#initialize sfs object
sfs = np.zeros((pop0Proj + 1, pop1Proj + 1))

#loop through variants
for variant in vcf:
    
    #get genotypes
    gts=variant.genotypes
    
    #get minor allele
    allele = 1
    
    #get minor allele counts
    minorAlleleCounts = getAlleleCounts(pop0Indices,pop1Indices,gts,allele)
    '''
    #get sample sizes
    sampSizePop0 = len([allele for gt in (gts[Idx][0:(len(gts[Idx]) - 1)] for Idx in pop0Indices) for allele in gt])
    sampSizePop1 = len([allele for gt in (gts[Idx][0:(len(gts[Idx]) - 1)] for Idx in pop1Indices) for allele in gt])
    '''
    #get sample sizes
    sampSizePop0 = len([allele for gt in (gts[Idx][0:(len(gts[Idx]) - 1)] for Idx in pop0Indices) for allele in gt if allele != -1])
    sampSizePop1 = len([allele for gt in (gts[Idx][0:(len(gts[Idx]) - 1)] for Idx in pop1Indices) for allele in gt if allele != -1])

    if(sampSizePop0 < pop0Proj or sampSizePop1 < pop1Proj):
        #print('skipping locus')
        continue
    
    #get projected probs for both populations
    projProbPop0 = projectCounts(minorAlleleCounts[0],sampSizePop0,pop0Proj)
    projProbPop1 = projectCounts(minorAlleleCounts[1],sampSizePop1,pop1Proj)
    
    #loop through projection values
    for i in range(len(projProbPop0)):
        for j in range(len(projProbPop1)):
            
            #calculate jiont probability
            jointProb = projProbPop0[i] * projProbPop1[j]
            
            #flip alleles depending on current projection bin
            if i + j > (pop0Proj + pop1Proj)/2:
                i_fold = pop0Proj - i
                j_fold = pop1Proj - j
            else:
                i_fold = i
                j_fold = j
            
            #add to sfs bin
            sfs[i_fold,j_fold] += jointProb
            
#Build SFS
sfsDf = pd.DataFrame(sfs)
#colnames and rownames
sfsDf.index = [('d' + args.pop1[3:4] + '_' + str(sampSize)) for sampSize in list(range(len(projProbPop0)))]
sfsDf.columns = [('d' + args.pop2[3:4] + '_' + str(sampSize)) for sampSize in list(range(len(projProbPop1)))]

#convert to tsv string
sfsTSV = sfsDf.to_csv(sep='\t')

#Save SFS
with open(args.outSFS, "w") as file:
    file.write("1 observation\n")
    file.write(sfsTSV)

#print statement
print("Saved " + args.outSFS)
    








