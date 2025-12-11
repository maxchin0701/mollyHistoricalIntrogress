#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 13:42:57 2025

@author: maxchin
"""

#import modules
from Bio import AlignIO
from Bio import SeqIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
import argparse
import pandas as pd
import sys

#function for getting LOH in alignment
#how do we communicate between alignments?
#pass on info from previous LOH
def findLOH(alignment,prevAlignLOH,prevAlignLOHAnc,progress,lastCheck,
            chromInterval,chromLength):
    
    #get current chromosome
    chrom=alignment[0].id.split(".")[1]
    #check if all samples are represented
    if(len(alignment) != 9):
        #index variable for current ungapped length
        ungapLength=0
        #loop through all positions
        for curPos in range(0,len(alignment[0].seq)):
            #progress meter
            #check if current position is gap
            if alignment[0].seq[curPos] != "-":
                #if current ungapped position is equal to position of next check
                if (progress + ungapLength) == (lastCheck + chromInterval) :
                    #print progress message
                    print(chrom, ": ", 
                          round((progress)/chromLength*100,1),
                          "% done",sep="")
                    #update check variable
                    lastCheck=lastCheck+chromInterval
                #increment ungapped length
                ungapLength+=1
        #update progress
        progress=progress+alignment[0].annotations["size"]
        
        #end block if previous alignment is unended
        if(prevAlignLOH!="NA"):
            #print("Current alignment incomplete, ending anc + modern block from last alignment")
            LOHList=prevAlignLOH
        else:
            LOHList="NA"
        
        #end block if previous alignment is unended
        if(prevAlignLOHAnc!="NA"):
            #print("Current alignment incomplete, ending anc block from last alignment")
            LOHListAnc=prevAlignLOHAnc
        else:
            LOHListAnc="NA"
        
        #return
        return(LOHList,LOHListAnc,progress,lastCheck,"NA","NA")
    
    #retreive samp idxs
    hapMexIdx=[k for k, curSamp in enumerate(alignment) if 
               curSamp.id.split(".")[0] == 'PForHapMex'][0]
    
    hapLatIdx=[k for k, curSamp in enumerate(alignment) if 
               curSamp.id.split(".")[0] == 'PForHaplat'][0]
    pMexIdx=[k for k, curSamp in enumerate(alignment) if 
             curSamp.id.split(".")[0] == 'PMex'][0]
    pLatIdx=[k for k, curSamp in enumerate(alignment) if 
             curSamp.id.split(".")[0] == 'PLat'][0]
    ancestralHapMexIdx=[k for k, curSamp in enumerate(alignment) if 
                        curSamp.id.split(".")[0] == 'AncPForHapMex'][0]
    ancestralHapLatIdx=[k for k, curSamp in enumerate(alignment) if 
                        curSamp.id.split(".")[0] == 'AncPForHapLat'][0]
    
    sampIdx=[hapMexIdx,hapLatIdx,pMexIdx,
             pLatIdx,ancestralHapMexIdx,ancestralHapLatIdx]
    
    #loop through and check for fixed differences, then LOH sites as well
    if(prevAlignLOH=="NA"):
        LOHList=[]
        prevLocLOH=False
        consecLOH=0
        ungapLength=0
    else:
        LOHList=prevAlignLOH
        prevLocLOH=True
        consecLOH=LOHList[0][3]-1
        ungapLength=0
        #print("Continuing Anc + modern LOH block from last locus")
        #print("Consec LOH:", consecLOH)
    
    #storage for LOH ancestral only
    if(prevAlignLOHAnc=="NA"):
        LOHListAnc=[]
        prevLocLOHAnc=False
        consecLOHAnc=0
    else:
        LOHListAnc=prevAlignLOHAnc
        prevLocLOHAnc=True
        consecLOHAnc=LOHListAnc[0][3]-1
        #print("Continuing Anc LOH block from last locus")
        #print("Consec Anc LOH:", consecLOHAnc)

    for curPos in range(0,len(alignment[hapMexIdx].seq)):
        #progress meter
        #check if current position is gap
        if alignment[0].seq[curPos] != "-":
            #if current ungapped position is equal to position of next check
            if (progress + ungapLength) == (lastCheck + chromInterval) :
                #print progress message
                print(chrom, ": ",
                      round((progress)/chromLength*100,1),
                      "% done",sep="")
                #update check variable
                lastCheck=lastCheck+chromInterval
            #increment ungapped length
            ungapLength+=1
            
        #check for fixed difference
        if(alignment[pMexIdx].seq[curPos].capitalize() != 
           alignment[pLatIdx].seq[curPos].capitalize() and
           "-" not in [([curSamp.seq[curPos] for curSamp in alignment])[idx] 
                       for idx in sampIdx]): 
            #print("fixedDiff")
            #check for LOH
            if(alignment[ancestralHapMexIdx].seq[curPos].capitalize() == 
              alignment[ancestralHapLatIdx].seq[curPos].capitalize() and 
              alignment[ancestralHapMexIdx].seq[curPos].capitalize() in 
              [alignment[pMexIdx].seq[curPos].capitalize(),
               alignment[pLatIdx].seq[curPos].capitalize()]):
                #print("Ancestral LOH at locus")
                #get ancestry at locus
                if(alignment[ancestralHapMexIdx].seq[curPos].capitalize() == 
                   alignment[pMexIdx].seq[curPos].capitalize()):
                    ancRetain="pmex"
                else:
                    ancRetain="plat"
                    
                #If previous fixed difference is not block
                if(prevLocLOHAnc==False or 
                   prevLocLOHAnc==True and 
                   consecLOHAnc < 1):
                    #update consecLOH
                    if(prevLocLOHAnc == True and 
                       LOHListAnc[len(LOHListAnc)-1][4] == 
                       ancRetain):
                        consecLOHAnc+=1
                    
                    #update prevLocLOH
                    prevLocLOHAnc=True
                    
                    #set nSNPs
                    nSNPs=1
                    
                    #append LOH
                    LOHListAnc.append([alignment[hapMexIdx].id.split(".")[1],
                                    (progress+ungapLength)-1,
                                    (progress+ungapLength)+1,
                                    nSNPs,
                                    ancRetain])
                #if at least 2 consecutive fixed diffs have been LOH
                elif(prevLocLOHAnc==True and 
                        consecLOHAnc >= 1 and
                        LOHListAnc[len(LOHListAnc)-1][4] == ancRetain):
                    #if initiating a new block
                    if(consecLOHAnc==1):
                        #print('starting Anc block')
                        #update ending position and nSNPs of new block
                        LOHListAnc[len(LOHListAnc)-2][2]=LOHListAnc[len(LOHListAnc)-1][2]
                        LOHListAnc[len(LOHListAnc)-2][3]+=1
                        #remove most recent LOH snp (now part of block)
                        del LOHListAnc[len(LOHListAnc)-1]
                    
                    #add to block
                    #print("adding to existing Anc block")
                    LOHListAnc[len(LOHListAnc)-1][2]=(progress+ungapLength)+1
                    LOHListAnc[len(LOHListAnc)-1][3]+=1
                    
                    #increment consecLOH
                    consecLOHAnc+=1
                    prevLocLOHAnc==True
                
                #print(*LOHList,sep='\n')
                
                if(alignment[hapMexIdx].seq[curPos].capitalize() == 
                  alignment[hapLatIdx].seq[curPos].capitalize() and 
                  alignment[hapMexIdx].seq[curPos].capitalize() == 
                    alignment[ancestralHapMexIdx].seq[curPos].capitalize() and
                  alignment[hapLatIdx].seq[curPos].capitalize() == 
                    alignment[ancestralHapLatIdx].seq[curPos].capitalize() and
                  alignment[hapMexIdx].seq[curPos].capitalize() in 
                  [alignment[pMexIdx].seq[curPos].capitalize(),
                   alignment[pLatIdx].seq[curPos].capitalize()]):
                    #print("Ancestral+modern LOH at locus")
                    #get ancestry at locus
                    if(alignment[hapMexIdx].seq[curPos].capitalize() == 
                       alignment[pMexIdx].seq[curPos].capitalize()):
                        ancRetain="pmex"
                    else:
                        ancRetain="plat"
                        
                    #If previous fixed difference is not 
                    if(prevLocLOH==False or 
                       prevLocLOH==True and 
                       consecLOH < 1):
                        #update consecLOH
                        if(prevLocLOH == True and 
                           LOHList[len(LOHList)-1][4] == 
                           ancRetain):
                            consecLOH+=1
                        
                        #update prevLocLOH
                        prevLocLOH=True
                        
                        #set nSNPs
                        nSNPs=1
                        
                        #append LOH
                        LOHList.append([alignment[hapMexIdx].id.split(".")[1],
                                        (progress+ungapLength)-1,
                                        (progress+ungapLength)+1,
                                        nSNPs,
                                        ancRetain])
                    #if at least 2 consecutive fixed diffs have been LOH
                    elif(prevLocLOH==True and 
                            consecLOH >= 1 and
                            LOHList[len(LOHList)-1][4] == ancRetain):
                        #if initiating a new block
                        if(consecLOH==1):
                            #print('starting Anc + modern block')
                            #update ending position and nSNPs of new block
                            LOHList[len(LOHList)-2][2]=LOHList[len(LOHList)-1][2]
                            LOHList[len(LOHList)-2][3]+=1
                            #remove most recent LOH snp (now part of block)
                            del LOHList[len(LOHList)-1]
                        
                        #add to block
                        LOHList[len(LOHList)-1][2]=(progress+ungapLength)+1
                        LOHList[len(LOHList)-1][3]+=1
                        
                        #increment consecLOH
                        consecLOH+=1
                        prevLocLOH==True

            #if current fixed diff is not LOH
            else:    
                #current fixed diff is not LOH
                prevLocLOH=False
                prevLocLOHAnc=False
                consecLOH=0
                consecLOHAnc=0
    
    #check if alignment ended with LOH
    if(prevLocLOH==True and consecLOH > 1):
        endLOH=[LOHList[len(LOHList)-1]]
        del LOHList[len(LOHList)-1]
    else:
        endLOH="NA"
        
    #check if alignment ended with ancestral only LOH
    if(prevLocLOHAnc==True and consecLOHAnc > 1):
        endLOHAnc=[LOHListAnc[len(LOHListAnc)-1]]
        del LOHListAnc[len(LOHListAnc)-1]
    else:
        endLOHAnc="NA"

    #update progress
    progress=progress+alignment[hapMexIdx].annotations["size"]
    return(LOHList,LOHListAnc,progress,lastCheck,endLOH,endLOHAnc)
    
#parse user inputs
parser = argparse.ArgumentParser(description="Parsing test MAF output")
parser.add_argument("--inMAF", action="store", dest="inMAF", type=str)
#parser.add_argument("--inFASTA", action="store", dest="inFASTA", type=str)
parser.add_argument("--inIndex", action="store", dest="inIndex", type=str)
parser.add_argument("--inChrs", action="store", dest="inChrs", type=str)
parser.add_argument("--outLOHRegions", action="store", dest="outLOHRegions", type=str)
parser.add_argument("--outLOHRegionsAnc", action="store", dest="outLOHRegionsAnc", type=str)
args=parser.parse_args()

#read in chromosomes
chrDf=pd.read_csv(args.inChrs, sep="\t")
chroms=chrDf.iloc[:,0]
#chroms=["chr2"]

#print # chrs
print("traversing",len(chroms),"chromosomes")

#list of all LOH tracts
LOHAll=[]
LOHAllAnc=[]

#loop through chromosomes
for chrom in chroms:
    #print chromosome
    print("traversing",chrom)
    #get length of chromosome and chromosome interval 
    chromLength=chrDf[chrDf["CHR"]==chrom]["END"].iloc[0]-1
    #print(chromLength)
    chromInterval=int(chromLength/1000)
    progress=0
    lastCheck=0
        
    #load in current chromosomes 
    curChromAlignments=list(AlignIO.parse(args.inMAF+"/"+chrom+".maf", "maf"))
    #print(len(curChromAlignments))
    #print(chromInterval)

    #list to store current chromosome LOH
    LOHCurChrom=[]
    LOHAncCurChrom=[]
    prevAlignLOH="NA"
    prevAlignLOHAnc="NA"
    #loop through alignments
    for curAlign in curChromAlignments:
        #print('new alignment')
        #parse alignment for LOH
        curAlignLOH=findLOH(curAlign,prevAlignLOH,prevAlignLOHAnc,
                            progress,lastCheck,
                            chromInterval,chromLength)
        #print(curAlignLOH[0])
        #print(*curAlignLOH[0],sep='\n')
        #print(len(curAlignLOH[0]))
        #extend current chromosome LOH
        if(curAlignLOH[0] != 'NA'):
            LOHCurChrom.extend(curAlignLOH[0])
            #print(len(LOHCurChrom))
        
        #extend current chromosome ancestral LOH
        if(curAlignLOH[1] != 'NA'):
            LOHAncCurChrom.extend(curAlignLOH[1])
            #print(len(LOHAncCurChrom))
        
        #update progress and checkpoints
        progress=curAlignLOH[2]
        lastCheck=curAlignLOH[3]
        prevAlignLOH=curAlignLOH[4]
        #print(prevAlignLOH)
        prevAlignLOHAnc=curAlignLOH[5]
        #print(prevAlignLOHAnc)
    
    print(chrom,"traversal complete")
    #print(*LOHCurChrom,sep='\n')
    #extend LOHAll
    LOHAll.extend(LOHCurChrom)
    LOHAllAnc.extend(LOHAncCurChrom)
    #print(len(LOHCurChrom))
    #print(len(LOHAncCurChrom))

#give titles  to df
LOHAllDf=pd.DataFrame(LOHAll,columns=["chr","start","end","nSNPs","anc"])
LOHAllDf["intSize"]=LOHAllDf["end"]-LOHAllDf["start"]
LOHAllAncDf=pd.DataFrame(LOHAllAnc,columns=["chr","start","end","nSNPs","anc"])
LOHAllAncDf["intSize"]=LOHAllAncDf["end"]-LOHAllAncDf["start"]

#print(LOHAllDf)
#print(LOHAllAncDf)

#write to tsv
LOHAllDf.to_csv(args.outLOHRegions,sep='\t',index=False)
LOHAllAncDf.to_csv(args.outLOHRegionsAnc,sep='\t',index=False)

#print(sum(LOHAllDf["intSize"])/len(LOHAllDf))
#print(sum(LOHAllAncDf["intSize"])/len(LOHAllAncDf))

        
    



