#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 16:49:37 2025

@author: maxchin
"""

import argparse
import pandas
pandas.options.mode.chained_assignment = None
import os
import statistics
import arviz as az
import numpy as np

#parse user inputs
parser = argparse.ArgumentParser(description="Summing CDS across LOH-GFF intersection")
parser.add_argument("--inIntersect", action="store", dest="inIntersect", type=str)
parser.add_argument("--inSimIntersect",action="store",dest="inSimIntersect",type=str)
args=parser.parse_args()

#define function for gettting longest transcripts
def getLongestTranscripts(inIntersect):
    #get list of genes
    genes=[]
    for i in range(0,len(inIntersect)):
        if inIntersect.iloc[i,5] == "gene": 
            curGene=inIntersect.iloc[i,11].split('=')[1]
            if curGene not in genes:
                genes.append(curGene)
    
    transcripts=[]
    dfMRNAs=inIntersect[inIntersect.iloc[:,5]=='mRNA']
    for gene in genes:
        curGeneMRNAs=dfMRNAs[dfMRNAs.iloc[:,11].str.split(
            ';',expand=True)[1].str.split('=',expand=True)[1] == gene]
        curGeneMRNAs['intSize']=curGeneMRNAs.iloc[:,7]-curGeneMRNAs.iloc[:,6]
        transcripts.append(curGeneMRNAs[curGeneMRNAs['intSize'] == max(curGeneMRNAs['intSize'])].iloc[
            0,11].split(';')[0].split('=')[1])
    return(transcripts)

def calcCDSOverlap(LOHStart,LOHEnd,CDSStart,CDSEnd): 
    if (CDSStart >= LOHStart and 
        CDSEnd <= LOHEnd): 
        return(CDSEnd-CDSStart)
    elif (CDSStart < LOHStart and 
          CDSEnd <= LOHEnd):
        return(CDSEnd-LOHStart)
    elif (CDSStart >= LOHStart and 
          CDSEnd > LOHEnd):
        return(LOHEnd-CDSStart)
    else :
        return(LOHEnd-LOHStart)


#parse through empirical df
dfIntersect = pandas.read_csv(args.inIntersect,sep='\t',header=None)
keepTranscriptsEmp=getLongestTranscripts(dfIntersect)

#store empirical lengths
sumCDSEmp=0

#sum empirical cds 
for i in range(0,len(dfIntersect)):
    print("processing empirical output")
    #isolate to just CDS rows
    if (dfIntersect.iloc[i,5] == "CDS" and 
        dfIntersect.iloc[i,11].split(';')[1].split('=')[1] in keepTranscriptsEmp): 
        sumCDSEmp+=calcCDSOverlap(dfIntersect.iloc[i,1], 
                                  dfIntersect.iloc[i,2],
                                  dfIntersect.iloc[i,6],
                                  dfIntersect.iloc[i,7])
        
#parse through sim outputs
sumCDSSimAll=[]
simIter=0
for simOut in os.scandir(args.inSimIntersect):
    simIter+=1
    print("analyzing simulation iteration ", simIter)
    #read in simulation output
    dfSimIntersect = pandas.read_csv(simOut,sep='\t',header=None)
    keepTranscriptsSim=getLongestTranscripts(dfSimIntersect)
    sumCDSSim=0
    for i in range(0,len(dfSimIntersect)):
        if (dfSimIntersect.iloc[i,5] == "CDS" and 
            dfSimIntersect.iloc[i,11].split(';')[1].split('=')[1] in keepTranscriptsSim): 
                sumCDSSim+=calcCDSOverlap(dfSimIntersect.iloc[i,1], 
                                          dfSimIntersect.iloc[i,2],
                                          dfSimIntersect.iloc[i,6],
                                          dfSimIntersect.iloc[i,7])
    #print(simOut,": ",sumCDSSim)
    #add to list
    sumCDSSimAll.append(sumCDSSim)

#calculate p-value
nGOE=0
for i in sumCDSSimAll:
    if i >= sumCDSEmp:
        nGOE+=1

print("Total empirical CDS: ",sumCDSEmp)
print("Average simulated total CDS: ",statistics.mean(sumCDSSimAll))
print("p-value: ",nGOE/1000)
print("HPD interval: ",az.hdi(np.array(sumCDSSimAll), hdi_prob=0.95))
print("Done")


    

