#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script saves gene ENSEMBL IDs and gene symbols for genes intersecting
empirical shared LOH tracts
"""

import argparse
import pandas as pd

#parse user inputs
parser = argparse.ArgumentParser(description="Extracting gene names from intersection")
parser.add_argument("--inIntersect", action="store", dest="inIntersect", type=str)
parser.add_argument("--outGenes",action="store",dest="outGenes",type=str)
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

#parse through empirical df
dfIntersect = pd.read_csv(args.inIntersect,sep='\t',header=None)
keepTranscriptsEmp=getLongestTranscripts(dfIntersect)

#store empirical lengths
geneN=[]
geneRefseq=[]
geneSymbols=[]

#sum empirical cds 
for i in range(0,len(dfIntersect)):
    #isolate to just mRNA rows
    if dfIntersect.iloc[i,5] == "mRNA" and dfIntersect.iloc[
            i,11].split(';')[0].split('=')[1] in keepTranscriptsEmp:
        #split to extract field identifiers
        idFields=[element.split('=')[0] for element in dfIntersect.iloc[
            i,11].split(';')]
        #skip if no accession
        if len(idFields) <= 2:
            continue
        #Append gene name
        geneN.append(dfIntersect.iloc[i,11].split(';')[0].split('=')[1])
        #append accession
        geneRefseq.append(dfIntersect.iloc[i,11].split(';')[2].split('=')[
            1].split('.')[1])        
        #check if symbol present
        if 'em_Preferred_name' in idFields:
            symbolIndex=[i for i, curField in enumerate(idFields) if 
                          curField == 'em_Preferred_name'][0]
            symbol=dfIntersect.iloc[i,11].split(';')[symbolIndex].split('=')[1]
            if symbol == "":
                symbol="NA"
        else:
            symbol="NA"
        #append symbol
        geneSymbols.append(symbol)
'''
#print number of genes
print(len(geneN))
print(geneRefseq)
'''
#save
out=open(args.outGenes,'w')
for i in range(len(geneN)):
    out.write(geneN[i]+'\t'+geneRefseq[i]+'\t'+geneSymbols[i]+'\n')
out.close()


            
            
            
            