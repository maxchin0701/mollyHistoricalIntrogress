# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 17:16:03 2025

@author: maxchin
"""

#import modules
from Bio import AlignIO
from Bio import SeqIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
import argparse

#parse user inputs
parser = argparse.ArgumentParser(description="Parsing test MAF output")
parser.add_argument("--inMAF", action="store", dest="inMAF", type=str)
parser.add_argument("--inFASTA", action="store", dest="inFASTA", type=str)
parser.add_argument("--inIndex", action="store", dest="inIndex", type=str)
parser.add_argument("--indexSamp", action="store", dest="indexSamp", type=str)
args=parser.parse_args()

testIndex = MafIO.MafIndex(args.inIndex,
                           args.inMAF,
                           "simHuman_chr6.simHuman.chr6")
"""
ungapLength=[{"species":"simHuman_chr6.simHuman.chr6","len":0},
             {"species":"mrh.mrhrefChr0","len":0},
             {"species":"mrhcd.mrhcdrefChr0","len":0},
             {"species":"mr.mrrefChr0","len":0},
             {"species":"simMouse_chr6.simMouse.chr6","len":0},
             {"species":"simRat_chr6.simRat.chr6","len":0},
             {"species":"cd.cdrefChr0","len":0},
             {"species":"simCow_chr6.simCow.chr6","len":0},
             {"species":"simDog_chr6.simDog.chr6","len":0}]
                        
                
totLength=[{"species":"simHuman_chr6.simHuman.chr6","len":0},
             {"species":"mrh.mrhrefChr0","len":0},
             {"species":"mrhcd.mrhcdrefChr0","len":0},
             {"species":"mr.mrrefChr0","len":0},
             {"species":"simMouse_chr6.simMouse.chr6","len":0},
             {"species":"simRat_chr6.simRat.chr6","len":0},
             {"species":"cd.cdrefChr0","len":0},
             {"species":"simCow_chr6.simCow.chr6","len":0},
             {"species":"simDog_chr6.simDog.chr6","len":0}]
"""
'''
numIncomplete=0
numAlignments=0

for align in AlignIO.parse(args.inMAF, "maf"):
    numAlignments+=1
    curLengths=[]
    if len(align) != 9:
        numIncomplete+=1
    """   
    for seqRec in align:
        ungapSeq=seqRec.seq.replace("-","")
        curSampIndex=[i for i, curSeqRec in enumerate(ungapLength) if 
                      curSeqRec["species"] == seqRec.id][0]
        #ungapLength[curSampIndex]["len"]+=1
        #ungapLength[curSampIndex]["len"]+=len(ungapSeq)
        curLengths.append(len(seqRec.seq))
        totLength[curSampIndex]["len"]+=len(seqRec.seq)
        #print(len(ungapSeq))
        #print(len(align[0].seq),"\n")
    curTotLengths=[]
    for i in range(len(totLength)):
        curTotLengths.append(totLength[i]["len"])
    
    #print(curTotLengths)
    print(curLengths)
    if (len(set(curLengths)) != 1 or
        len(curLengths) != 9):
        print(align)
    """
print(numAlignments)
print(numIncomplete/numAlignments)
#print(ungapLength)
#print(totLength)
'''
testIntersect=testIndex.search([100567],[100593])
print(testIntersect)

for align in testIntersect:
    print("new alignment\n\n")
    print(align,"\n")

"""    
    #get record length and start positions
    for seqRec in align:
        print(seqRec.id)
        print(seqRec.seq)
        print(seqRec.annotations["start"])
        print(seqRec.annotations["start"] + seqRec.annotations["size"],"\n")

for record in SeqIO.parse(args.inFASTA, "fasta"):
    print(record.seq[100458:100734])
"""




